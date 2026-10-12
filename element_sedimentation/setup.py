#!/usr/bin/env python3
"""
Install the element_sedimentation extension into a MESA tree.

The extension consists of five Fortran sources that are copied into
$MESA_DIR/star/private and a handful of MESA files that are patched by anchor
insertion. Each patch snippet in ./code consists of blocks

    ! anchor = '<text of a line in the MESA file>'[; replace | ; prepend | ; all]
    <lines to insert>

The lines are inserted after the (unique) MESA line containing <text>, or
instead of it (replace), before it (prepend), or after every occurrence (all).
Pristine MESA files are backed up to ./backup together with the MESA version
they came from; uninstall restores them.

    python3 setup.py <command> [-f NAME] [-y] [--no-build] [--force] [-v] [--mesa-dir DIR]

Commands:
    install            patch MESA, copy the sources, full MESA rebuild (./clean; ./install)
    soft_install       patch MESA, copy the sources, rebuild star_data and star only
    uninstall          restore MESA, remove the sources, rebuild star_data and star
    soft_uninstall     same as uninstall (kept for compatibility)
    reinstall          uninstall + install (full rebuild)
    soft_reinstall     uninstall + soft_install
    single_install     patch/copy one file (-f NAME) and rebuild its module
    single_uninstall   restore/remove one file (-f NAME) and rebuild its module
    single_reinstall   restore + patch one file (-f NAME) and rebuild its module
    status             which MESA files are patched, which sources are present, versions
    check              dry run: parse every snippet and verify the anchors (MESA is not modified)
    clean              remove the backup directory

Importing this module has no side effects (MESA_DIR is only read when needed).
Requires Python >= 3.9, standard library only.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
CODE_DIR = SCRIPT_DIR / "code"
BACKUP_DIR = SCRIPT_DIR / "backup"
VERSION_FILE = BACKUP_DIR / "MESA_VERSION"

# a MESA file is pristine if it does not contain this substring
MARKER = "element_sedimentation"

# the new sources (copied, in dependency order) and the patched MESA files
NEW_SOURCES = (
    "advection_diffusion.f90",
    "sedimentation_phase_diagram.f90",
    "sedimentation_miscibility.f90",
    "element_sedimentation_legacy.f90",
    "element_sedimentation.f90",
)

INSTALL_INFO: dict[str, dict] = {
    # star/make
    "makefile_base": {"mesa_path": "star/make", "replace": False},
    # star/defaults
    "controls.defaults": {"mesa_path": "star/defaults", "replace": False},
    # star/private (patched)
    "alloc.f90": {"mesa_path": "star/private", "replace": False},
    "ctrls_io.f90": {"mesa_path": "star/private", "replace": False},
    "evolve.f90": {"mesa_path": "star/private", "replace": False},
    "hydro_energy.f90": {"mesa_path": "star/private", "replace": False},
    "init.f90": {"mesa_path": "star/private", "replace": False},
    # star_data
    "star_controls.inc": {"mesa_path": "star_data/private", "replace": False},
    "star_data_step_work.inc": {"mesa_path": "star_data/public", "replace": False},
}
for _name in NEW_SOURCES:
    INSTALL_INFO[_name] = {"mesa_path": "star/private", "replace": True}

# MESA modules to rebuild after patching (order matters: star depends on star_data;
# ./clean is required because the .inc changes alter the star_info type)
BUILD_DIRS = ("star_data", "star")

MODE_APPEND = "append"
MODE_REPLACE = "replace"
MODE_PREPEND = "prepend"
MODE_APPEND_ALL = "all"

ANCHOR_RE = re.compile(r"^! anchor = '(.*)'\s*(?:;\s*(replace|prepend|all))?\s*$")

VERBOSE = False


class SetupError(Exception):
    """A problem that should be reported to the user without a traceback."""


def log(msg: str) -> None:
    if VERBOSE:
        print(msg)


def hr() -> None:
    print("-" * 80)


# ---------------------------------------------------------------------------
# MESA tree
# ---------------------------------------------------------------------------

_MESA_DIR_OVERRIDE: Path | None = None


def mesa_dir() -> Path:
    """Path of the MESA installation (--mesa-dir or $MESA_DIR)."""
    if _MESA_DIR_OVERRIDE is not None:
        path = _MESA_DIR_OVERRIDE
    else:
        value = os.environ.get("MESA_DIR")
        if not value:
            raise SetupError(
                "the environment variable MESA_DIR is not set "
                "(point it to your MESA installation or use --mesa-dir)"
            )
        path = Path(value)
    if not path.is_dir():
        raise SetupError(f"MESA_DIR = {path} is not a directory")
    return path


def mesa_version() -> str:
    """Contents of $MESA_DIR/data/version_number ('unknown' if absent)."""
    path = mesa_dir() / "data" / "version_number"
    if not path.is_file():
        return "unknown"
    return path.read_text().strip()


def backup_version() -> str | None:
    """MESA version the backups were taken from, None if no version file."""
    if not VERSION_FILE.is_file():
        return None
    return VERSION_FILE.read_text().strip()


def is_patched_text(text: str) -> bool:
    return MARKER in text


# ---------------------------------------------------------------------------
# snippets
# ---------------------------------------------------------------------------


@dataclass
class Snippet:
    anchor: str
    mode: str
    lines: list[str] = field(default_factory=list)


def parse_snippet_text(text: str, name: str) -> list[Snippet]:
    """Split a snippet file into (anchor, mode, lines) blocks."""
    snippets: list[Snippet] = []
    for lineno, line in enumerate(text.splitlines(keepends=True), start=1):
        if line.startswith("! anchor ="):
            match = ANCHOR_RE.match(line.rstrip("\r\n"))
            if match is None:
                raise SetupError(f"{name}:{lineno}: malformed anchor line: {line.rstrip()}")
            anchor, mode = match.group(1), match.group(2) or MODE_APPEND
            if not anchor:
                raise SetupError(f"{name}:{lineno}: empty anchor")
            snippets.append(Snippet(anchor, mode))
        elif snippets:
            snippets[-1].lines.append(line)
        elif line.strip():
            raise SetupError(f"{name}:{lineno}: content before the first anchor")
    if not snippets:
        raise SetupError(f"no anchor found in {name}")
    for sn in snippets:
        if not sn.lines and sn.mode != MODE_REPLACE:
            raise SetupError(f"{name}: anchor {sn.anchor!r} has no content")
        if sn.lines and not sn.lines[-1].endswith("\n"):
            sn.lines[-1] += "\n"
    return snippets


def parse_snippet_file(path: Path) -> list[Snippet]:
    if not path.is_file():
        raise SetupError(f"snippet file {path} not found")
    return parse_snippet_text(path.read_text(), path.name)


def find_anchor(lines: list[str], sn: Snippet, name: str) -> list[int]:
    """Indices of the lines containing the anchor; raises unless the count is valid."""
    hits = [i for i, line in enumerate(lines) if sn.anchor in line]
    if not hits:
        raise SetupError(f"anchor {sn.anchor!r} not found in {name}")
    if len(hits) > 1 and sn.mode != MODE_APPEND_ALL:
        raise SetupError(
            f"anchor {sn.anchor!r} found {len(hits)} times in {name} "
            "(use '; all' to patch every occurrence)"
        )
    return hits


def apply_snippets(lines: list[str], snippets: list[Snippet], name: str) -> list[str]:
    """Return a patched copy of lines (lines of a pristine MESA file)."""
    lines = list(lines)
    for sn in snippets:
        for i in reversed(find_anchor(lines, sn, name)):
            if sn.mode in (MODE_APPEND, MODE_APPEND_ALL):
                start, end = i + 1, i + 1
            elif sn.mode == MODE_REPLACE:
                start, end = i, i + 1
            elif sn.mode == MODE_PREPEND:
                start, end = i, i
            else:  # pragma: no cover - the regex only admits the modes above
                raise SetupError(f"unknown insert mode {sn.mode!r}")
            if sn.mode == MODE_PREPEND:
                present = lines[max(i - len(sn.lines), 0) : i] == sn.lines
            else:
                present = bool(sn.lines) and lines[start : start + len(sn.lines)] == sn.lines
            if present:
                log(f"   content of anchor {sn.anchor!r} already present in {name}")
                continue
            lines[start:end] = sn.lines
    return lines


# ---------------------------------------------------------------------------
# installation files
# ---------------------------------------------------------------------------


class InstallationFile:
    """One file of the extension: a snippet patching a MESA file, or a new source."""

    def __init__(self, name: str) -> None:
        if name not in INSTALL_INFO:
            raise SetupError(f"unknown file {name!r}; known files: {', '.join(INSTALL_INFO)}")
        info = INSTALL_INFO[name]
        self.name = name
        self.replace: bool = info["replace"]
        self.code_path = CODE_DIR / name
        self.mesa_subdir = Path(info["mesa_path"])
        self.mesa_file = mesa_dir() / self.mesa_subdir / name
        self.module_dir = mesa_dir() / self.mesa_subdir.parts[0]  # star or star_data
        self.backup_path = BACKUP_DIR / (name + ".bak")
        if not self.code_path.is_file():
            raise SetupError(f"{self.code_path} not found")
        if not self.replace and not self.mesa_file.is_file():
            raise SetupError(f"{self.mesa_file} not found")

    # --- state -------------------------------------------------------------

    def is_patched(self) -> bool:
        """True if the MESA file contains the extension (for new sources: exists)."""
        if self.replace:
            return self.mesa_file.is_file()
        return is_patched_text(self.mesa_file.read_text())

    def has_backup(self) -> bool:
        return self.backup_path.is_file()

    @property
    def backup_version_path(self) -> Path:
        return BACKUP_DIR / (self.name + ".bak.version")

    def backup_version(self) -> str | None:
        """MESA version this file's backup was taken from (falls back to the global record)."""
        if self.backup_version_path.is_file():
            return self.backup_version_path.read_text().strip()
        return backup_version()

    def backup_text(self) -> str:
        """Text of the backup; refuses a backup that is not pristine."""
        if not self.has_backup():
            raise SetupError(f"no backup found for {self.name} in {BACKUP_DIR}")
        text = self.backup_path.read_text()
        if is_patched_text(text):
            raise SetupError(
                f"the backup {self.backup_path} is not pristine (it contains '{MARKER}');\n"
                "   restore the MESA file by hand (e.g. from a fresh MESA download), then delete the backup"
            )
        return text

    def is_up_to_date(self) -> bool | None:
        """True if the MESA file is exactly what installing the current code/ version gives,
        False if it differs, None if that cannot be decided (no usable pristine backup)."""
        if self.replace:
            return self.mesa_file.is_file() and self.mesa_file.read_bytes() == self.code_path.read_bytes()
        if not self.is_patched():
            return False
        try:
            lines = self.pristine_text().splitlines(keepends=True)
        except SetupError:
            return None
        try:
            expected = "".join(apply_snippets(lines, parse_snippet_file(self.code_path), self.name))
        except SetupError:
            return False
        return self.mesa_file.read_text() == expected

    def pristine_text(self) -> str:
        """Text of the pristine MESA file (from MESA itself, or from the backup)."""
        if self.replace:
            raise SetupError(f"{self.name} is a new source file; it has no pristine version")
        text = self.mesa_file.read_text()
        if not is_patched_text(text):
            return text
        if self.has_backup():
            return self.backup_text()
        raise SetupError(f"{self.mesa_file} is patched and no backup exists in {BACKUP_DIR}")

    # --- actions -----------------------------------------------------------

    def backup(self) -> None:
        """(Re)take the backup whenever the MESA file is pristine."""
        if self.replace:
            return
        if not self.is_patched():
            BACKUP_DIR.mkdir(parents=True, exist_ok=True)
            shutil.copy2(self.mesa_file, self.backup_path)
            # per-file record only: the global backup/MESA_VERSION is kept as a read-only
            # fallback for backups made by older versions of this script, so retaking one
            # backup must not change the recorded origin of the others
            self.backup_version_path.write_text(mesa_version() + "\n")
            log(f"   backed up {self.name} to {self.backup_path}")
        elif not self.has_backup():
            raise SetupError(
                f"{self.mesa_file} is already patched but no backup exists in {BACKUP_DIR};\n"
                "   restore the pristine file by hand (e.g. from a fresh MESA download) before installing"
            )

    def check_backup_version(self, force: bool) -> None:
        current, backed = mesa_version(), self.backup_version()
        if backed is None:
            msg = f"the backup of {self.name} has no MESA version record"
        elif backed != current:
            msg = f"the backup of {self.name} is from MESA {backed} but MESA_DIR is {current}"
        else:
            return
        if force:
            print(f"WARNING: {msg}; restoring anyway (--force)")
        else:
            raise SetupError(msg + " (use --force to restore anyway)")

    def restore(self, force: bool = False) -> None:
        """Copy the backup over the (patched) MESA file."""
        if not self.has_backup():
            raise SetupError(f"no backup found for {self.name} in {BACKUP_DIR}")
        self.check_backup_version(force)
        text = self.backup_text()
        # write (do not copy2): the restored file must be newer than its object file,
        # otherwise a build without ./clean keeps the patched object
        self.mesa_file.write_text(text)
        log(f"   restored {self.name}")

    def install(self, force: bool = False) -> None:
        """Copy a new source, or patch the MESA file (from its pristine version)."""
        if self.replace:
            shutil.copy(self.code_path, self.mesa_file)  # current mtime, so that make recompiles it
            log(f"   copied {self.name} into {self.mesa_file.parent}")
            return
        self.backup()
        if self.is_patched():
            # patch the pristine version again: this makes install idempotent and
            # picks up changes of the snippet
            self.restore(force)
            # the MESA file is pristine now: retake the backup, which records its MESA version
            # (after a --force restore of a backup without a record this adopts it)
            self.backup()
        snippets = parse_snippet_file(self.code_path)
        lines = self.mesa_file.read_text().splitlines(keepends=True)
        lines = apply_snippets(lines, snippets, self.name)
        self.mesa_file.write_text("".join(lines))
        log(f"   patched {self.name}")

    def uninstall(self, force: bool = False) -> None:
        """Remove a new source, or restore the MESA file from the backup."""
        if self.replace:
            if self.mesa_file.is_file():
                self.mesa_file.unlink()
                log(f"   removed {self.mesa_file}")
            return
        if not self.is_patched():
            log(f"   {self.name} is pristine, nothing to restore")
            return
        self.restore(force)

    def check(self) -> list[str]:
        """Dry run: problems found while applying the snippet to the pristine file."""
        problems: list[str] = []
        if self.replace:
            return problems
        try:
            snippets = parse_snippet_file(self.code_path)
            lines = self.pristine_text().splitlines(keepends=True)
            for sn in snippets:
                hits = find_anchor(lines, sn, self.name)
                log(f"   {self.name}: anchor {sn.anchor!r} ({sn.mode}) found {len(hits)}x")
            apply_snippets(lines, snippets, self.name)
        except SetupError as exc:
            problems.append(str(exc))
        return problems


def all_files() -> list[InstallationFile]:
    return [InstallationFile(name) for name in INSTALL_INFO]


# ---------------------------------------------------------------------------
# building
# ---------------------------------------------------------------------------


def run(cmd: list[str], cwd: Path) -> None:
    """Run one of MESA's build scripts (./clean, ./mk, ./export, ./install) in cwd.

    MESA's scripts are not all executable with a shebang line (e.g. star_data/clean
    is a plain list of shell commands), so they are always run through bash.
    """
    print(f"$ (cd {cwd} && {' '.join(cmd)})")
    script = cwd / cmd[0]
    if not script.exists():
        raise SetupError(f"cannot run {cmd[0]} in {cwd}: no such file")
    env = dict(os.environ)
    env["MESA_DIR"] = str(mesa_dir())  # --mesa-dir must also apply to MESA's build scripts
    try:
        subprocess.run(["bash", str(script)] + cmd[1:], cwd=cwd, check=True, env=env)
    except OSError as exc:
        raise SetupError(f"cannot run {cmd[0]} in {cwd}: {exc}") from exc
    except subprocess.CalledProcessError as exc:
        raise SetupError(f"{' '.join(cmd)} in {cwd} failed with exit status {exc.returncode}") from exc


def build_modules(dirs=BUILD_DIRS, clean: bool = True) -> None:
    """./clean; ./mk; ./export in the given MESA module directories."""
    for name in dirs:
        cwd = mesa_dir() / name
        if clean:
            run(["./clean"], cwd)
        run(["./mk"], cwd)
        run(["./export"], cwd)


def full_rebuild() -> None:
    """./clean; ./install in MESA_DIR (rebuilds every MESA module)."""
    run(["./clean"], mesa_dir())
    run(["./install"], mesa_dir())


def confirm(question: str, assume_yes: bool) -> bool:
    if assume_yes:
        return True
    if not sys.stdin.isatty():
        raise SetupError(f"{question} -- no terminal to confirm; use -y/--yes")
    return input(f"{question} [y/N]: ").strip().lower() in ("y", "yes")


def remind() -> None:
    hr()
    print("Don't forget to run ./clean and ./mk in your work directory.")
    hr()


def check_backup_consistency(force: bool, files: list[InstallationFile] | None = None) -> None:
    """Before anything is modified: every patched MESA file must have a pristine backup taken
    from the current MESA version (version mismatches can be overridden with --force).
    Checking up front keeps a failing install/uninstall from leaving MESA half restored."""
    current = mesa_version()
    problems, stale = [], []
    for f in files if files is not None else all_files():
        if f.replace or not f.is_patched():
            continue
        try:
            f.backup_text()
        except SetupError as exc:
            problems.append(str(exc))
            continue
        backed = f.backup_version()
        if backed is None or backed != current:
            stale.append(f"{f.name} (backup from MESA {backed or 'unknown'})")
    if problems:
        raise SetupError("cannot restore the pristine MESA files:\n   " + "\n   ".join(problems))
    if stale:
        msg = (f"MESA_DIR ({current}) is patched but these backups are not recorded as taken from it: "
               f"{', '.join(stale)}")
        if not force:
            raise SetupError(
                msg + ".\n   If the backups are from this MESA (e.g. made by an older version of this "
                "script, which did not record versions), run 'python3 setup.py soft_install --force' once: "
                "it restores and retakes them with a version record. Otherwise restore the files by hand.")
        print(f"WARNING: {msg}; restoring anyway (--force)")


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------


def check_snippets_apply(files: list[InstallationFile] | None = None) -> None:
    """Dry run of every snippet against the pristine MESA files; raise before anything is
    modified if an anchor is missing (e.g. after a MESA upgrade)."""
    problems = [msg for f in (files if files is not None else all_files()) for msg in f.check()]
    if problems:
        raise SetupError("the patch snippets do not apply to this MESA; nothing was modified:\n   "
                         + "\n   ".join(problems))


def do_install_files(force: bool) -> None:
    check_backup_consistency(force)
    check_snippets_apply()
    for f in all_files():
        print(f"installing {f.name}")
        f.install(force)


def do_uninstall_files(force: bool) -> None:
    check_backup_consistency(force)
    for f in all_files():
        print(f"uninstalling {f.name}")
        f.uninstall(force)


def cmd_install(args) -> None:
    if not args.no_build and not confirm(
        f"Patch MESA at {mesa_dir()} and run a full rebuild (./clean; ./install)?", args.yes
    ):
        raise SetupError("aborted")
    do_install_files(args.force)
    if not args.no_build:
        full_rebuild()
    hr()
    print("element_sedimentation was successfully installed.")
    remind()


def cmd_soft_install(args) -> None:
    do_install_files(args.force)
    if not args.no_build:
        build_modules()
    hr()
    print("element_sedimentation was successfully installed (star_data and star rebuilt).")
    remind()


def cmd_uninstall(args) -> None:
    do_uninstall_files(args.force)
    if not args.no_build:
        build_modules()
    hr()
    print("element_sedimentation was successfully uninstalled (star_data and star rebuilt).")
    remind()


def cmd_reinstall(args) -> None:
    if not args.no_build and not confirm(
        f"Re-patch MESA at {mesa_dir()} and run a full rebuild (./clean; ./install)?", args.yes
    ):
        raise SetupError("aborted")
    check_backup_consistency(args.force)
    check_snippets_apply()  # a reinstall that cannot succeed must not uninstall first
    do_uninstall_files(args.force)
    do_install_files(args.force)
    if not args.no_build:
        full_rebuild()
    hr()
    print("element_sedimentation was successfully reinstalled.")
    remind()


def cmd_soft_reinstall(args) -> None:
    check_backup_consistency(args.force)
    check_snippets_apply()  # a reinstall that cannot succeed must not uninstall first
    do_uninstall_files(args.force)
    do_install_files(args.force)
    if not args.no_build:
        build_modules()
    hr()
    print("element_sedimentation was successfully reinstalled (star_data and star rebuilt).")
    remind()


def single_file(args) -> InstallationFile:
    name = args.filename
    if not name:
        if args.yes or not sys.stdin.isatty():
            raise SetupError("single_* commands need -f/--filename NAME")
        name = input("Enter the file name: ").strip()
    return InstallationFile(name)


def build_for_single(f: InstallationFile, args) -> None:
    if args.no_build:
        return
    if f.module_dir.name == "star_data":
        # star_info changed: star must be rebuilt from scratch as well
        build_modules(("star_data", "star"))
    else:
        build_modules(("star",), clean=False)


def cmd_single_install(args) -> None:
    f = single_file(args)
    check_backup_consistency(args.force, [f])
    check_snippets_apply([f])
    f.install(args.force)
    build_for_single(f, args)
    hr()
    print(f"{f.name} was successfully installed.")
    remind()


def cmd_single_uninstall(args) -> None:
    f = single_file(args)
    check_backup_consistency(args.force, [f])
    f.uninstall(args.force)
    build_for_single(f, args)
    hr()
    print(f"{f.name} was successfully uninstalled.")
    remind()


def cmd_single_reinstall(args) -> None:
    f = single_file(args)
    check_backup_consistency(args.force, [f])
    check_snippets_apply([f])
    f.uninstall(args.force)
    f.install(args.force)
    build_for_single(f, args)
    hr()
    print(f"{f.name} was successfully reinstalled.")
    remind()


def cmd_status(args) -> int:
    print(f"MESA_DIR:       {mesa_dir()}")
    print(f"MESA version:   {mesa_version()}")
    current_version = mesa_version()
    if BACKUP_DIR.is_dir():
        print(f"backup dir:     {BACKUP_DIR}")
    else:
        print(f"backup dir:     {BACKUP_DIR} (absent)")
    print()
    print("patched MESA files:")
    n_patched = 0
    n_stale = 0
    n_unknown = 0
    version_mismatch = []
    files = all_files()
    for f in files:
        if f.replace:
            continue
        patched = f.is_patched()
        n_patched += patched
        state = "pristine"
        if patched:
            current = f.is_up_to_date()
            if current is None:
                n_unknown += 1
                state = "patched (no pristine backup)"
            elif current:
                state = "patched (up to date)"
            else:
                n_stale += 1
                state = "patched (STALE)"
        if f.has_backup():
            backed = f.backup_version()
            backup_info = f"backup from {backed or 'unknown MESA'}"
            if backed != current_version:
                version_mismatch.append(f.name)
        else:
            backup_info = "no backup"
        print(f"   {str(f.mesa_subdir / f.name):40s} {state:28s} {backup_info}")
    print()
    print(f"new sources in {mesa_dir() / 'star' / 'private'}:")
    n_present = 0
    for f in files:
        if not f.replace:
            continue
        present = f.is_patched()
        n_present += present
        state = "missing"
        if present:
            current = f.is_up_to_date()
            n_stale += not current
            state = "present (up to date)" if current else "present (STALE)"
        print(f"   {f.name:40s} {state}")
    n_patch_files = sum(1 for f in files if not f.replace)
    print()
    if n_patched == 0 and n_present == 0:
        print("status: not installed")
    elif n_patched == n_patch_files and n_present == len(NEW_SOURCES):
        if n_unknown > 0:
            print(f"status: installed, but {n_unknown} patched file(s) have no pristine backup in "
                  f"{BACKUP_DIR}; they cannot be verified, updated or uninstalled. Restore them from a "
                  "fresh MESA download, then run install.")
        elif n_stale == 0:
            print("status: installed and up to date with code/")
        else:
            print(f"status: installed but {n_stale} file(s) differ from code/ -- run "
                  "'python3 setup.py soft_install' to update MESA")
    else:
        print("status: PARTIALLY installed -- run 'python3 setup.py soft_install' or 'uninstall'")
    if version_mismatch:
        print(f"WARNING: backups not recorded as taken from MESA {current_version}: {', '.join(version_mismatch)}")
    return 0


def cmd_check(args) -> int:
    """Dry run of the anchor parser against the pristine MESA files; returns the exit status."""
    print(f"MESA_DIR: {mesa_dir()} (version {mesa_version()})")
    problems: list[str] = []
    missing = [name for name in INSTALL_INFO if not (CODE_DIR / name).is_file()]
    for name in missing:
        problems.append(f"{CODE_DIR / name} is missing")
    for name in INSTALL_INFO:
        if name in missing:
            continue
        try:
            f = InstallationFile(name)
        except SetupError as exc:
            problems.append(str(exc))
            continue
        if f.replace:
            print(f"   {name:28s} new source, {f.code_path.stat().st_size} bytes")
            continue
        file_problems = f.check()
        snippets = parse_snippet_file(f.code_path) if not file_problems else []
        source = "MESA file" if not is_patched_text(f.mesa_file.read_text()) else "backup"
        state = "OK" if not file_problems else "FAILED"
        print(f"   {name:28s} {len(snippets):2d} anchor(s) checked against the {source}: {state}")
        problems.extend(file_problems)
    problems.extend(check_control_consistency())
    print()
    if problems:
        print("check FAILED:")
        for p in problems:
            print(f"   - {p}")
        return 1
    print("check passed: every anchor was found, every snippet applies cleanly, "
          "and the controls are consistent")
    return 0


def check_control_consistency() -> list[str]:
    """Every control declared in star_controls.inc must be initialized in controls.defaults
    and appear in the namelist, the 's% x = x' and the 'x = s% x' blocks of ctrls_io.f90
    (MESA includes controls.defaults at compile time, so a missing initializer leaves the
    control undefined without any compiler diagnostic)."""
    problems: list[str] = []
    try:
        decl = (CODE_DIR / "star_controls.inc").read_text()
        defaults = (CODE_DIR / "controls.defaults").read_text()
        io = (CODE_DIR / "ctrls_io.f90").read_text()
    except OSError as exc:
        return [f"cannot read the control snippets: {exc}"]
    names: list[str] = []
    joined = ""
    for raw in decl.splitlines():
        if raw.startswith("! anchor ="):
            continue
        line = raw.split("!", 1)[0].rstrip()  # no strings in declarations, so '!' starts a comment
        if not line.strip():
            continue
        joined += line.strip()
        if joined.endswith("&"):  # continuation line
            joined = joined[:-1] + " "
            continue
        stmt, joined = joined, ""
        if "::" not in stmt:
            problems.append(f"star_controls.inc: cannot parse declaration {stmt!r}")
            continue
        for item in stmt.split("::", 1)[1].split(","):
            match = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)", item)
            if match is None:
                problems.append(f"star_controls.inc: cannot parse declaration {stmt!r}")
                break
            names.append(match.group(1))
    if not names:
        return problems + ["no control declarations found in star_controls.inc"]
    for name in names:
        if not re.search(rf"^\s*{name}\s*=", defaults, flags=re.M):
            problems.append(f"control {name} is declared in star_controls.inc but not initialized in controls.defaults")
        if not re.search(rf"^\s*{name},\s*&", io, flags=re.M):
            problems.append(f"control {name} is missing from the namelist block of ctrls_io.f90")
        if not re.search(rf"^\s*s%\s*{name}\s*=\s*{name}\s*$", io, flags=re.M):
            problems.append(f"control {name} is missing from the 's% x = x' block of ctrls_io.f90")
        if not re.search(rf"^\s*{name}\s*=\s*s%\s*{name}\s*$", io, flags=re.M):
            problems.append(f"control {name} is missing from the 'x = s% x' block of ctrls_io.f90")
    for name in re.findall(r"^\s*(element_sedimentation_\w+|do_element_sedimentation)\s*=", defaults, flags=re.M):
        if name not in names:
            problems.append(f"control {name} is initialized in controls.defaults but not declared in star_controls.inc")
    print(f"   {len(names)} controls checked for consistency across star_controls.inc, controls.defaults and ctrls_io.f90")
    return problems


def cmd_clean(args) -> None:
    if not BACKUP_DIR.is_dir():
        print(f"{BACKUP_DIR} does not exist, nothing to clean")
        return
    patched = [f.name for f in all_files() if not f.replace and f.is_patched()]
    if patched and not args.force:
        raise SetupError(
            f"MESA is still patched ({', '.join(patched)}); uninstall first, or use --force to "
            "delete the backups anyway"
        )
    if not confirm(f"Delete {BACKUP_DIR}?", args.yes):
        raise SetupError("aborted")
    shutil.rmtree(BACKUP_DIR)
    print(f"removed {BACKUP_DIR}")


COMMANDS = {
    "install": cmd_install,
    "uninstall": cmd_uninstall,
    "reinstall": cmd_reinstall,
    "clean": cmd_clean,
    "soft_install": cmd_soft_install,
    "soft_uninstall": cmd_uninstall,
    "soft_reinstall": cmd_soft_reinstall,
    "single_install": cmd_single_install,
    "single_uninstall": cmd_single_uninstall,
    "single_reinstall": cmd_single_reinstall,
    "status": cmd_status,
    "check": cmd_check,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Install or uninstall the element_sedimentation extension in $MESA_DIR.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Commands:", 1)[1],
    )
    parser.add_argument("command", choices=list(COMMANDS), help="the command to execute")
    parser.add_argument("-f", "--filename", help="file name for the single_* commands")
    parser.add_argument("-y", "--yes", action="store_true", help="do not ask for confirmation")
    parser.add_argument("--no-build", action="store_true", help="patch/restore only, do not build")
    parser.add_argument("--force", action="store_true",
                        help="ignore the MESA version check of the backups; allow clean while patched")
    parser.add_argument("-v", "--verbose", action="store_true", help="more output")
    parser.add_argument("--mesa-dir", help="MESA installation to use instead of $MESA_DIR")
    return parser


def main(argv: list[str] | None = None) -> int:
    global VERBOSE, _MESA_DIR_OVERRIDE
    sys.stdout.reconfigure(line_buffering=True)  # keep progress and error messages in order
    args = build_parser().parse_args(argv)
    VERBOSE = args.verbose
    if args.mesa_dir:
        _MESA_DIR_OVERRIDE = Path(args.mesa_dir).expanduser().resolve()
    try:
        if not CODE_DIR.is_dir():
            raise SetupError(f"{CODE_DIR} not found")
        result = COMMANDS[args.command](args)
    except SetupError as exc:
        print(f"setup.py {args.command}: error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
    return int(result or 0)


if __name__ == "__main__":
    sys.exit(main())
