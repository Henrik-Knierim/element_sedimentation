#!/usr/bin/env python3
"""
Install the element_sedimentation extension into a MESA tree.

The extension consists of Fortran sources that are copied into
$MESA_DIR/star/private and a handful of MESA files that are patched by anchor
insertion. Each patch snippet in ./code consists of blocks

    ! anchor = '<text of a line in the MESA file>'[; replace | ; prepend | ; all]
    <lines to insert>

The lines are inserted after the (unique) MESA line containing <text>, or
instead of it (replace), before it (prepend), or after every occurrence (all).

The script records exactly what it inserted (./backup/<file>.applied.json) and
removes exactly that again on uninstall, so several extensions that use this
installer (element_sedimentation, gentle_mixing) can patch the same MESA files
in any order and be updated or removed independently. The MESA file as found
before the first install is also kept in ./backup together with the MESA
version it came from; it is only used with --force, when nothing else works.

    python3 setup.py <command> [-f NAME] [-y] [--no-build] [--force] [-v] [--mesa-dir DIR]

Commands:
    install            patch MESA, copy the sources, full MESA rebuild (./clean; ./install)
    soft_install       patch MESA, copy the sources, rebuild star_data and star only
    uninstall          unpatch MESA, remove the sources, rebuild star_data and star
    soft_uninstall     same as uninstall (kept for compatibility)
    reinstall          uninstall + install (full rebuild)
    soft_reinstall     uninstall + soft_install
    single_install     patch/copy one file (-f NAME) and rebuild its module
    single_uninstall   unpatch/remove one file (-f NAME) and rebuild its module
    single_reinstall   unpatch + patch one file (-f NAME) and rebuild its module
    status             which MESA files are patched, which sources are present, versions
    check              dry run: parse every snippet and verify the anchors (MESA is not modified)
    clean              remove the backup directory (backups and records)

Importing this module has no side effects (MESA_DIR is only read when needed).
Requires Python >= 3.9, standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ===========================================================================
# configuration of this extension: the only part of this script that differs
# between the extensions that share it (element_sedimentation, gentle_mixing)
# ===========================================================================

EXTENSION = "element_sedimentation"

# a MESA file is patched by this extension if it contains this substring
# (every block of every snippet must contain it)
MARKER = "element_sedimentation"

# the new sources (copied into star/private, in dependency order)
NEW_SOURCES = (
    "advection_diffusion.f90",
    "sedimentation_phase_diagram.f90",
    "sedimentation_miscibility.f90",
    "sedimentation_energy.f90",
    "sedimentation_overturn.f90",
    "element_sedimentation.f90",
)

# sources of earlier versions of the extension that no longer exist; install and uninstall
# remove them from star/private so that an updated checkout leaves no stray file behind
RETIRED_SOURCES = (
    "element_sedimentation_legacy.f90",  # the old instant_rain / bottom_up_rain, removed 2026-10-14
)

# the patched MESA files: snippet name -> directory in MESA_DIR
PATCHED_FILES = {
    "makefile_base": "star/make",
    "controls.defaults": "star/defaults",
    "alloc.f90": "star/private",
    "ctrls_io.f90": "star/private",
    "evolve.f90": "star/private",
    "hydro_energy.f90": "star/private",
    "init.f90": "star/private",
    "star_controls.inc": "star_data/private",
    "star_data_step_work.inc": "star_data/public",
}

# MESA modules to rebuild after patching (order matters: star depends on star_data;
# ./clean is required because the .inc changes alter the star_info type)
BUILD_DIRS = ("star_data", "star")

# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CODE_DIR = SCRIPT_DIR / "code"
BACKUP_DIR = SCRIPT_DIR / "backup"
VERSION_FILE = BACKUP_DIR / "MESA_VERSION"

INSTALL_INFO: dict[str, dict] = {name: {"mesa_path": sub, "replace": False} for name, sub in PATCHED_FILES.items()}
for _name in NEW_SOURCES:
    INSTALL_INFO[_name] = {"mesa_path": "star/private", "replace": True}

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


@dataclass
class AppliedBlock:
    """One block as it was inserted into a MESA file (what uninstall removes again)."""
    anchor: str
    mode: str
    lines: list[str]
    replaced: list[str] = field(default_factory=list)  # the original line, mode replace only
    count: int = 1  # occurrences (mode all: one per anchor hit)

    def key(self) -> tuple:
        return (self.anchor, self.mode, tuple(self.lines))

    def to_json(self) -> dict:
        return {"anchor": self.anchor, "mode": self.mode, "lines": self.lines,
                "replaced": self.replaced, "count": self.count}

    @staticmethod
    def from_json(d: dict) -> "AppliedBlock":
        return AppliedBlock(d["anchor"], d["mode"], list(d["lines"]), list(d.get("replaced", [])),
                            int(d.get("count", 1)))


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
        if not sn.lines:
            raise SetupError(f"{name}: anchor {sn.anchor!r} has no content")
        if not sn.lines[-1].endswith("\n"):
            sn.lines[-1] += "\n"
        if not any(MARKER in line for line in sn.lines):
            raise SetupError(
                f"{name}: the block of anchor {sn.anchor!r} does not contain the marker {MARKER!r} "
                "(add it in a comment; status detects patched files by it)")
    return snippets


def parse_snippet_file(path: Path) -> list[Snippet]:
    if not path.is_file():
        raise SetupError(f"snippet file {path} not found")
    return parse_snippet_text(path.read_text(), path.name)


def find_anchor(lines: list[str], anchor: str, mode: str, name: str, block: list[str] | None = None) -> list[int]:
    """Indices of the lines containing the anchor (copies of the anchor text inside an already
    inserted copy of `block` do not count); raises unless the count is valid."""
    inside: set[int] = set()
    if block:
        for j in _find_block(lines, block):
            inside.update(range(j, j + len(block)))
    hits = [i for i, line in enumerate(lines) if anchor in line and i not in inside]
    if not hits:
        raise SetupError(f"anchor {anchor!r} not found in {name}")
    if len(hits) > 1 and mode != MODE_APPEND_ALL:
        raise SetupError(
            f"anchor {anchor!r} found {len(hits)} times in {name} "
            "(use '; all' to patch every occurrence)"
        )
    return hits


def apply_snippets(lines: list[str], snippets: list[Snippet], name: str) -> tuple[list[str], list[AppliedBlock]]:
    """Return a patched copy of lines and the record of the blocks inserted."""
    lines = list(lines)
    applied: list[AppliedBlock] = []
    for sn in snippets:
        hits = find_anchor(lines, sn.anchor, sn.mode, name, sn.lines)
        block = AppliedBlock(sn.anchor, sn.mode, list(sn.lines), count=len(hits))
        for i in reversed(hits):
            if sn.mode in (MODE_APPEND, MODE_APPEND_ALL):
                start, end = i + 1, i + 1
            elif sn.mode == MODE_REPLACE:
                start, end = i, i + 1
                block.replaced = [lines[i]]
            elif sn.mode == MODE_PREPEND:
                start, end = i, i
            else:  # pragma: no cover - the regex only admits the modes above
                raise SetupError(f"unknown insert mode {sn.mode!r}")
            if sn.mode == MODE_PREPEND:
                present = lines[max(i - len(sn.lines), 0) : i] == sn.lines
            else:
                present = lines[start : start + len(sn.lines)] == sn.lines
            if present:
                print(f"   NOTE: the content of anchor {sn.anchor!r} is already present in {name}; "
                      "it is recorded as this extension's and will be removed on uninstall")
                if sn.mode == MODE_REPLACE:
                    block.replaced = []
                continue
            lines[start:end] = sn.lines
        applied.append(block)
    return lines, applied


def _find_block(lines: list[str], block_lines: list[str], lo: int = 0, hi: int | None = None) -> list[int]:
    """Start indices j in [lo, hi) with lines[j:j+n] == block_lines."""
    n = len(block_lines)
    hi = len(lines) if hi is None else hi
    return [j for j in range(lo, max(lo, hi - n + 1)) if lines[j : j + n] == block_lines]


def remove_blocks(lines: list[str], blocks: list[AppliedBlock], name: str) -> list[str]:
    """Return a copy of lines without the recorded blocks (the reverse of apply_snippets).

    Blocks are located by their content, not by position: another extension may have
    inserted its own lines after the same anchor, before or after ours."""
    lines = list(lines)
    for block in reversed(blocks):
        n = len(block.lines)
        if block.mode == MODE_REPLACE:
            hits = _find_block(lines, block.lines)
            if len(hits) != 1:
                raise SetupError(
                    f"{name}: the block recorded for anchor {block.anchor!r} was found {len(hits)} times "
                    "(expected once); was the file edited by hand? Restore it from the backup with --force")
            lines[hits[0] : hits[0] + n] = block.replaced
            continue
        # anchor lines, not counting copies of the anchor text inside the inserted block itself
        # (a block may legitimately repeat the line it was anchored to)
        hits = _find_block(lines, block.lines)
        anchors = [i for i, line in enumerate(lines)
                   if block.anchor in line and not any(j <= i < j + n for j in hits)]
        if len(anchors) != block.count:
            raise SetupError(
                f"{name}: anchor {block.anchor!r} of a recorded block was found {len(anchors)} times "
                f"(expected {block.count}); was the file edited by hand? Restore it from the backup with --force")
        for i in reversed(anchors):
            hits = _find_block(lines, block.lines)
            if block.mode == MODE_PREPEND:
                candidates = [j for j in hits if j + n <= i]
                j = max(candidates) if candidates else None
            else:
                candidates = [j for j in hits if j > i]
                j = min(candidates) if candidates else None
            if j is None:
                raise SetupError(
                    f"{name}: the block recorded for anchor {block.anchor!r} was not found "
                    f"{'before' if block.mode == MODE_PREPEND else 'after'} the anchor; was the file edited by hand? "
                    "Restore it from the backup with --force")
            del lines[j : j + n]
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
        self.record_path = BACKUP_DIR / (name + ".applied.json")
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
        """Text of the backup; refuses a backup that contains this extension."""
        if not self.has_backup():
            raise SetupError(f"no backup found for {self.name} in {BACKUP_DIR}")
        text = self.backup_path.read_text()
        if is_patched_text(text):
            raise SetupError(
                f"the backup {self.backup_path} is not pristine (it contains '{MARKER}');\n"
                "   restore the MESA file by hand (e.g. from a fresh MESA download), then delete the backup"
            )
        return text

    def record(self) -> list[AppliedBlock] | None:
        """The blocks this script inserted, None if there is no record."""
        if not self.record_path.is_file():
            return None
        try:
            data = json.loads(self.record_path.read_text())
            return [AppliedBlock.from_json(d) for d in data["blocks"]]
        except (ValueError, KeyError, TypeError) as exc:
            raise SetupError(f"cannot read {self.record_path}: {exc}") from exc

    def write_record(self, blocks: list[AppliedBlock]) -> None:
        BACKUP_DIR.mkdir(parents=True, exist_ok=True)
        self.record_path.write_text(json.dumps(
            {"extension": EXTENSION, "file": self.name, "mesa_version": mesa_version(),
             "blocks": [b.to_json() for b in blocks]}, indent=1) + "\n")

    def synthesized_record(self) -> list[AppliedBlock] | None:
        """For a file patched by an older version of this script (no record): the blocks,
        reconstructed from the pristine backup and the current snippets, if the file is
        exactly backup + snippets. None otherwise."""
        if not self.has_backup():
            return None
        try:
            lines = self.backup_text().splitlines(keepends=True)
            patched, blocks = apply_snippets(lines, parse_snippet_file(self.code_path), self.name)
        except SetupError:
            return None
        if "".join(patched) != self.mesa_file.read_text():
            return None
        return blocks

    def blocks(self) -> list[AppliedBlock] | None:
        """Recorded or reconstructed blocks of this extension in the MESA file; None if the
        file is pristine or if its content cannot be accounted for."""
        rec = self.record()
        if rec is not None:
            return rec
        if not self.is_patched():
            return None
        return self.synthesized_record()

    def base_lines(self, force: bool = False) -> list[str]:
        """The MESA file without this extension's blocks (other extensions' patches stay)."""
        text = self.mesa_file.read_text()
        lines = text.splitlines(keepends=True)
        blocks = self.blocks()
        if blocks is not None:
            return remove_blocks(lines, blocks, self.name)
        if not is_patched_text(text):
            return lines
        # patched, but neither a record nor a reconstruction: the backup is the last resort
        if force and self.has_backup():
            self.check_backup_version(force)
            print(f"WARNING: {self.mesa_file} contains '{MARKER}' but no record of what was inserted; "
                  f"restoring the backup {self.backup_path} (--force). Patches of other extensions "
                  "in this file are lost; reinstall them afterwards.")
            return self.backup_text().splitlines(keepends=True)
        raise SetupError(
            f"{self.mesa_file} contains '{MARKER}' but there is no record in {BACKUP_DIR} of what this "
            "script inserted, and the file is not exactly backup + snippets.\n"
            "   Either restore the file by hand, or run again with --force to restore the backup "
            "(which also removes the patches of any other extension in this file).")

    def is_up_to_date(self) -> bool | None:
        """True if the MESA file contains exactly the blocks of the current code/ version,
        False if it differs, None if that cannot be decided."""
        if self.replace:
            return self.mesa_file.is_file() and self.mesa_file.read_bytes() == self.code_path.read_bytes()
        if not self.is_patched():
            return False
        blocks = self.blocks()
        if blocks is None:
            return None
        try:
            snippets = parse_snippet_file(self.code_path)
            remove_blocks(self.mesa_file.read_text().splitlines(keepends=True), blocks, self.name)
        except SetupError:
            return False
        return [b.key() for b in blocks] == [(sn.anchor, sn.mode, tuple(sn.lines)) for sn in snippets]

    def state(self) -> str:
        """One-line description of the patch state for status."""
        if not self.is_patched():
            return "pristine"
        if self.record() is None:
            if self.synthesized_record() is not None:
                return "patched by an older setup.py (migrates on next install)"
            return "patched, NO RECORD (cannot be unpatched; see --force)"
        current = self.is_up_to_date()
        if current:
            return "patched (up to date)"
        try:
            remove_blocks(self.mesa_file.read_text().splitlines(keepends=True), self.record(), self.name)
        except SetupError:
            return "patched (recorded blocks NOT FOUND: edited by hand?)"
        return "patched (STALE)"

    # --- actions -----------------------------------------------------------

    def backup(self) -> None:
        """Keep the MESA file as found before this extension first touched it."""
        if self.replace or self.has_backup():
            return
        if self.is_patched():
            log(f"   {self.name} already contains the extension and has no backup; none taken")
            return
        BACKUP_DIR.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self.mesa_file, self.backup_path)
        # per-file record only: the global backup/MESA_VERSION is kept as a read-only
        # fallback for backups made by older versions of this script
        self.backup_version_path.write_text(mesa_version() + "\n")
        log(f"   backed up {self.name} to {self.backup_path}")

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
        """Overwrite the MESA file with its backup (the file as found before the first install)."""
        self.check_backup_version(force)
        text = self.backup_text()
        if self.mesa_file.read_text() == text:
            return
        # a plain copy would keep the old mtime and make MESA's make skip the file
        self.mesa_file.write_text(text)
        log(f"   restored {self.name}")

    def install(self, force: bool = False) -> None:
        """Copy a new source, or (re)insert the snippets into the MESA file."""
        if self.replace:
            shutil.copy(self.code_path, self.mesa_file)  # current mtime, so that make recompiles it
            log(f"   copied {self.name} into {self.mesa_file.parent}")
            return
        self.backup()
        base = self.base_lines(force)  # the previous blocks of this extension removed
        patched, blocks = apply_snippets(base, parse_snippet_file(self.code_path), self.name)
        text = "".join(patched)
        if text != self.mesa_file.read_text():
            self.mesa_file.write_text(text)
            log(f"   patched {self.name}")
        else:
            log(f"   {self.name} is unchanged")
        self.write_record(blocks)

    def uninstall(self, force: bool = False) -> None:
        """Remove a new source, or remove this extension's blocks from the MESA file."""
        if self.replace:
            if self.mesa_file.is_file():
                self.mesa_file.unlink()
                log(f"   removed {self.mesa_file}")
            return
        if not self.is_patched() and self.record() is None:
            log(f"   {self.name} is pristine, nothing to do")
            return
        text = "".join(self.base_lines(force))
        if text != self.mesa_file.read_text():
            self.mesa_file.write_text(text)
            log(f"   unpatched {self.name}")
        if self.record_path.is_file():
            self.record_path.unlink()

    def check(self) -> list[str]:
        """Dry run: problems found while applying the snippets to the MESA file without this
        extension's blocks (other extensions' patches stay in place for the check)."""
        problems: list[str] = []
        if self.replace:
            return problems
        try:
            snippets = parse_snippet_file(self.code_path)
            lines = self.base_lines()
            for sn in snippets:
                hits = find_anchor(lines, sn.anchor, sn.mode, self.name, sn.lines)
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

    MESA's ./install also runs the test program of every module; some of them (colors)
    crash with the default OpenMP thread stack on macOS, and a failed ./install leaves
    the tree without libraries. OMP_STACKSIZE is therefore set for the build scripts
    unless the environment already defines it.
    """
    print(f"$ (cd {cwd} && {' '.join(cmd)})")
    script = cwd / cmd[0]
    if not script.exists():
        raise SetupError(f"cannot run {cmd[0]} in {cwd}: no such file")
    env = dict(os.environ)
    env["MESA_DIR"] = str(mesa_dir())  # --mesa-dir must also apply to MESA's build scripts
    env.setdefault("OMP_STACKSIZE", "512m")
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


def check_unpatchable(force: bool, files: list[InstallationFile] | None = None) -> None:
    """Before anything is modified: every patched MESA file must be unpatchable (a record,
    a reconstruction, or with --force a backup). Checking up front keeps a failing
    install/uninstall from leaving MESA half done."""
    problems = []
    for f in files if files is not None else all_files():
        if f.replace:
            continue
        try:
            f.base_lines(force)
        except SetupError as exc:
            problems.append(str(exc))
    if problems:
        raise SetupError("cannot unpatch the MESA files; nothing was modified:\n   " + "\n   ".join(problems))


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------


def check_snippets_apply(files: list[InstallationFile] | None = None) -> None:
    """Dry run of every snippet against the MESA files; raise before anything is
    modified if an anchor is missing (e.g. after a MESA upgrade)."""
    problems = [msg for f in (files if files is not None else all_files()) for msg in f.check()]
    if problems:
        raise SetupError("the patch snippets do not apply to this MESA; nothing was modified:\n   "
                         + "\n   ".join(problems))


def retired_sources_present() -> list[Path]:
    private = mesa_dir() / "star" / "private"
    return [private / name for name in RETIRED_SOURCES if (private / name).is_file()]


def remove_retired_sources() -> None:
    """Remove sources of earlier versions of the extension from star/private."""
    for path in retired_sources_present():
        path.unlink()
        print(f"removed the retired source {path}")


def do_install_files(force: bool) -> None:
    check_unpatchable(force)
    check_snippets_apply()
    for f in all_files():
        print(f"installing {f.name}")
        f.install(force)
    remove_retired_sources()


def do_uninstall_files(force: bool) -> None:
    check_unpatchable(force)
    for f in all_files():
        print(f"uninstalling {f.name}")
        f.uninstall(force)
    remove_retired_sources()


def cmd_install(args) -> None:
    if not args.no_build and not confirm(
        f"Patch MESA at {mesa_dir()} and run a full rebuild (./clean; ./install)?", args.yes
    ):
        raise SetupError("aborted")
    do_install_files(args.force)
    if not args.no_build:
        full_rebuild()
    hr()
    print(f"{EXTENSION} was successfully installed.")
    remind()


def cmd_soft_install(args) -> None:
    do_install_files(args.force)
    if not args.no_build:
        build_modules()
    hr()
    print(f"{EXTENSION} was successfully installed ({' and '.join(BUILD_DIRS)} rebuilt).")
    remind()


def cmd_uninstall(args) -> None:
    do_uninstall_files(args.force)
    if not args.no_build:
        build_modules()
    hr()
    print(f"{EXTENSION} was successfully uninstalled ({' and '.join(BUILD_DIRS)} rebuilt).")
    remind()


def cmd_reinstall(args) -> None:
    if not args.no_build and not confirm(
        f"Re-patch MESA at {mesa_dir()} and run a full rebuild (./clean; ./install)?", args.yes
    ):
        raise SetupError("aborted")
    check_unpatchable(args.force)
    check_snippets_apply()  # a reinstall that cannot succeed must not uninstall first
    do_uninstall_files(args.force)
    do_install_files(args.force)
    if not args.no_build:
        full_rebuild()
    hr()
    print(f"{EXTENSION} was successfully reinstalled.")
    remind()


def cmd_soft_reinstall(args) -> None:
    check_unpatchable(args.force)
    check_snippets_apply()  # a reinstall that cannot succeed must not uninstall first
    do_uninstall_files(args.force)
    do_install_files(args.force)
    if not args.no_build:
        build_modules()
    hr()
    print(f"{EXTENSION} was successfully reinstalled ({' and '.join(BUILD_DIRS)} rebuilt).")
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
    check_unpatchable(args.force, [f])
    check_snippets_apply([f])
    f.install(args.force)
    build_for_single(f, args)
    hr()
    print(f"{f.name} was successfully installed.")
    remind()


def cmd_single_uninstall(args) -> None:
    f = single_file(args)
    check_unpatchable(args.force, [f])
    f.uninstall(args.force)
    build_for_single(f, args)
    hr()
    print(f"{f.name} was successfully uninstalled.")
    remind()


def cmd_single_reinstall(args) -> None:
    f = single_file(args)
    check_unpatchable(args.force, [f])
    check_snippets_apply([f])
    f.uninstall(args.force)
    f.install(args.force)
    build_for_single(f, args)
    hr()
    print(f"{f.name} was successfully reinstalled.")
    remind()


def cmd_status(args) -> int:
    print(f"extension:      {EXTENSION}")
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
    n_migrate = 0
    version_mismatch = []
    files = all_files()
    for f in files:
        if f.replace:
            continue
        patched = f.is_patched()
        n_patched += patched
        state = f.state()
        if patched:
            if "NO RECORD" in state:
                n_unknown += 1
            elif "older setup.py" in state:
                n_migrate += 1
            elif state != "patched (up to date)":
                n_stale += 1
        if f.has_backup():
            backed = f.backup_version()
            backup_info = f"backup from {backed or 'unknown MESA'}"
            if backed != current_version:
                version_mismatch.append(f.name)
        else:
            backup_info = "no backup"
        print(f"   {str(f.mesa_subdir / f.name):40s} {state:52s} {backup_info}")
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
    for path in retired_sources_present():
        n_stale += 1
        print(f"   {path.name:40s} present (RETIRED, no longer part of the extension)")
    n_patch_files = sum(1 for f in files if not f.replace)
    print()
    if n_patched == 0 and n_present == 0:
        print("status: not installed")
    elif n_patched == n_patch_files and n_present == len(NEW_SOURCES):
        if n_unknown > 0:
            print(f"status: installed, but {n_unknown} patched file(s) have no record of what was inserted "
                  "and no usable backup; they cannot be updated or uninstalled automatically. "
                  "Restore them from a fresh MESA download, then run install.")
        elif n_stale == 0 and n_migrate > 0:
            print(f"status: installed and up to date with code/, but {n_migrate} file(s) were patched by an "
                  "older setup.py -- run 'python3 setup.py soft_install --no-build' once to record the blocks")
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
    """Dry run of the anchor parser against the MESA files; returns the exit status."""
    print(f"extension: {EXTENSION}")
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
        source = "MESA file" if not f.is_patched() else "MESA file without this extension's blocks"
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
        # namelist entries: one per line ('name, &') or several per line ('a, b, c, &')
        if not re.search(rf"(?<![\w%]){name}\s*,", io):
            problems.append(f"control {name} is missing from the namelist block of ctrls_io.f90")
        if not re.search(rf"^\s*s%\s*{name}\s*=\s*{name}\s*$", io, flags=re.M):
            problems.append(f"control {name} is missing from the 's% x = x' block of ctrls_io.f90")
        if not re.search(rf"^\s*{name}\s*=\s*s%\s*{name}\s*$", io, flags=re.M):
            problems.append(f"control {name} is missing from the 'x = s% x' block of ctrls_io.f90")
    for name in re.findall(r"^\s*([A-Za-z_]\w*)\s*=", defaults, flags=re.M):
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
            "delete the backups and records anyway"
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
        description=f"Install or uninstall the {EXTENSION} extension in $MESA_DIR.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Commands:", 1)[1],
    )
    parser.add_argument("command", choices=list(COMMANDS), help="the command to execute")
    parser.add_argument("-f", "--filename", help="file name for the single_* commands")
    parser.add_argument("-y", "--yes", action="store_true", help="do not ask for confirmation")
    parser.add_argument("--no-build", action="store_true", help="patch/unpatch only, do not build")
    parser.add_argument("--force", action="store_true",
                        help="restore the backup of a patched file that has no record (removes other "
                             "extensions' patches of it); ignore the MESA version check of the backups; "
                             "allow clean while patched")
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
