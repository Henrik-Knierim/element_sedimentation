#!/usr/bin/env python3
"""Run an example work directory with several element_sedimentation methods and compare.

    python3 compare_methods.py <example_dir> [options]

For every method a copy of <example_dir> is made in <out>/<method>/ (LOGS, png and
photos are not copied), ./mk is run, an inlist_run with the overrides is written and
read last by the top-level inlist, and ./star is run with the output in out.txt. The
history files are then compared at the end and at equal ages, and a table with the
wall-clock time per model, retries, termination, energy error, helium of the envelope
and centre, and the number of steps in which the module acted is printed and written
to <out>/results.json. The exit code is nonzero when a run did not reach the requested
age (unless --allow-failure).

Options
  --methods a,b,c     methods to run (default: all three)
  --max-age YEARS     stop age (default: taken from the example's inlist)
  --max-model N       max_model_number (default: as in the example)
  --load-model FILE   start from this saved model instead of the example's one (a path,
                      or the name of a model file in the example directory)
  --threads N         OMP_NUM_THREADS per run (default 4); runs are sequential
  --parallel          run all methods at once (needs enough cores)
  --out DIR           output directory (default: compare_methods_out next to the example)
  --extra 'name = value'   extra &controls line (repeatable)
  --allow-failure     exit 0 even if a run stops early
  --no-build          skip ./mk (the copies must have a star executable)
"""
import json
import os
import shutil
import subprocess
import sys
import time

METHODS = ['advection_diffusion_rain', 'instant_rain', 'bottom_up_rain']
SKIP = {'LOGS', 'png', 'photos', 'star', 'compare_methods_out', '.DS_Store'}


def read_history(path):
    with open(path) as f:
        lines = f.readlines()
    names = lines[5].split()
    rows = [[float(v) for v in l.split()] for l in lines[6:] if l.strip()]
    cols = {n: i for i, n in enumerate(names)}
    # keep the last entry of every model number (retries/backups)
    keep = {}
    for i, r in enumerate(rows):
        keep[int(r[cols['model_number']])] = i
    rows = [rows[i] for i in sorted(keep.values())]
    return cols, rows


def parse_args(argv):
    opts = {'methods': METHODS, 'max_age': None, 'max_model': None, 'load_model': None,
            'threads': 4, 'parallel': False, 'out': None, 'extra': [], 'allow_failure': False,
            'build': True}
    example = None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--methods':
            opts['methods'] = argv[i + 1].split(','); i += 2
        elif a == '--max-age':
            opts['max_age'] = float(argv[i + 1]); i += 2
        elif a == '--max-model':
            opts['max_model'] = int(argv[i + 1]); i += 2
        elif a == '--load-model':
            opts['load_model'] = argv[i + 1]; i += 2
        elif a == '--threads':
            opts['threads'] = int(argv[i + 1]); i += 2
        elif a == '--parallel':
            opts['parallel'] = True; i += 1
        elif a == '--out':
            opts['out'] = os.path.abspath(argv[i + 1]); i += 2
        elif a == '--extra':
            opts['extra'].append(argv[i + 1]); i += 2
        elif a == '--allow-failure':
            opts['allow_failure'] = True; i += 1
        elif a == '--no-build':
            opts['build'] = False; i += 1
        elif a.startswith('-'):
            sys.exit('unknown option ' + a)
        else:
            example = os.path.abspath(a); i += 1
    if example is None:
        sys.exit(__doc__)
    if opts['out'] is None:
        opts['out'] = os.path.join(os.path.dirname(example), 'compare_methods_out')
    if opts['load_model']:
        # a path relative to the current directory, or a file name in the example directory
        path = opts['load_model']
        if not os.path.isfile(path) and os.path.isfile(os.path.join(example, path)):
            path = os.path.join(example, path)
        if not os.path.isfile(path):
            sys.exit('--load-model: %s is neither a file nor a file in %s' % (opts['load_model'], example))
        opts['load_model'] = os.path.abspath(path)
    return example, opts


def copy_example(example, dest):
    if os.path.exists(dest):
        shutil.rmtree(dest)
    os.makedirs(dest)
    for name in os.listdir(example):
        if name in SKIP:
            continue
        src = os.path.join(example, name)
        dst = os.path.join(dest, name)
        if os.path.isdir(src):
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns('*.o', '*.mod', '*.smod'))
        else:
            shutil.copy2(src, dst)


def write_inlist_run(dest, method, opts):
    lines = ['&star_job', '   pgstar_flag = .false.']
    if opts['load_model']:
        lines += ['   load_saved_model = .true.',
                  "   load_model_filename = '%s'" % opts['load_model'],
                  '   set_initial_age = .false.',
                  '   set_initial_model_number = .false.',
                  '   set_initial_cumulative_energy_error = .true.']
    lines += ['/', '&controls',
              '   history_interval = 1',
              '   photo_interval = 100000',
              '   write_profiles_flag = .false.',
              '   element_sedimentation_report = .true.',
              '   do_element_sedimentation = .true.',
              "   element_sedimentation_method = '%s'" % method]
    if opts['max_age'] is not None:
        lines.append('   max_age = %.6e' % opts['max_age'])
    if opts['max_model'] is not None:
        lines.append('   max_model_number = %d' % opts['max_model'])
    for e in opts['extra']:
        lines.append('   ' + e)
    lines += ['/', '&pgstar', '/', '']
    with open(os.path.join(dest, 'inlist_run'), 'w') as f:
        f.write('\n'.join(lines))
    # read it last from the top-level inlist (after whatever the example reads)
    p = os.path.join(dest, 'inlist')
    s = open(p).read()
    for group, key in (('star_job', 'star_job'), ('controls', 'controls'), ('pgstar', 'pgstar')):
        n = 1
        while ('read_extra_%s_inlist(%d)' % (key, n)) in s:
            n += 1
        add = ('      read_extra_%s_inlist(%d) = .true.\n      extra_%s_inlist_name(%d) = \'inlist_run\'\n'
               % (key, n, key, n))
        marker = '&%s' % group
        i = s.index(marker) + len(marker)
        # insert after the group header line
        j = s.index('\n', i) + 1
        s = s[:j] + add + s[j:]
    open(p, 'w').write(s)


def run_one(example, method, opts):
    dest = os.path.join(opts['out'], method)
    copy_example(example, dest)
    write_inlist_run(dest, method, opts)
    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = str(opts['threads'])
    env.setdefault('OMP_STACKSIZE', '512m')
    if opts['build']:
        with open(os.path.join(dest, 'build.log'), 'w') as log:
            r = subprocess.run(['./mk'], cwd=dest, stdout=log, stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0 or not os.path.exists(os.path.join(dest, 'star')):
            return {'method': method, 'status': 'build failed', 'dir': dest}
    t0 = time.time()
    with open(os.path.join(dest, 'out.txt'), 'w') as log:
        subprocess.run(['./star'], cwd=dest, stdout=log, stderr=subprocess.STDOUT, env=env)
    return summarize(dest, method, time.time() - t0)


def summarize(dest, method, wall):
    res = {'method': method, 'dir': dest, 'wall_s': wall}
    out = open(os.path.join(dest, 'out.txt'), errors='replace').read()
    term = [l for l in out.splitlines() if l.startswith('termination code:')]
    res['termination'] = term[-1].split(':', 1)[1].strip() if term else 'none (crash?)'
    for l in reversed(out.splitlines()):
        if 'rel_cumulative_energy_error' in l:
            try:
                res['rel_cum_E_err'] = float(l.split()[-1].replace('D', 'E'))
            except ValueError:
                pass
            break
    lines = out.splitlines()
    res['module_steps'] = sum(1 for l in lines if l.startswith(' element_sedimentation: model')
                              and 'nothing to do' not in l and 'wall time' not in l)
    times = []
    for l in lines:
        if l.startswith(' element_sedimentation: model') and 'wall time of the module' in l:
            try:
                times.append(float(l.split()[-2]))
            except ValueError:
                pass
    res['module_wall_s_mean'] = sum(times) / len(times) if times else None
    res['module_failures'] = sum(1 for l in out.splitlines() if 'element_sedimentation' in l and ('failed' in l or 'error' in l.lower()) and 'report' not in l)
    hist = os.path.join(dest, 'LOGS', 'history.data')
    if os.path.exists(hist):
        cols, rows = read_history(hist)
        last = rows[-1]
        res['models'] = int(last[cols['model_number']])
        res['age_end'] = last[cols['star_age']]
        res['Teff_end'] = last[cols['Teff']]
        res['retries'] = int(last[cols['num_retries']]) if 'num_retries' in cols else None
        ycol = 'Y_env' if 'Y_env' in cols else ('surface_he4' if 'surface_he4' in cols else None)
        res['Y_env_end'] = last[cols[ycol]] if ycol else None
        res['Y_center_end'] = last[cols['center_he4']] if 'center_he4' in cols else None
        res['models_run'] = len(rows)
        res['wall_per_model_s'] = wall / max(len(rows), 1)
        res['history'] = {'age': [r[cols['star_age']] for r in rows],
                          'Teff': [r[cols['Teff']] for r in rows],
                          'Y_env': [r[cols[ycol]] for r in rows] if ycol else None,
                          'Y_center': [r[cols['center_he4']] for r in rows] if 'center_he4' in cols else None}
    else:
        res['models'] = 0
    return res


def interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if xs[i] >= x:
            f = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + f * (ys[i] - ys[i - 1])
    return ys[-1]


def main():
    example, opts = parse_args(sys.argv[1:])
    os.makedirs(opts['out'], exist_ok=True)
    print('example: %s\nmethods: %s\nout: %s' % (example, ', '.join(opts['methods']), opts['out']))
    results = []
    if opts['parallel']:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=len(opts['methods'])) as ex:
            results = list(ex.map(lambda m: run_one(example, m, opts), opts['methods']))
    else:
        for m in opts['methods']:
            print('running', m, '...', flush=True)
            results.append(run_one(example, m, opts))
    # table
    print()
    print('%-26s %7s %9s %9s %9s %8s %7s %9s %8s %8s %8s  %s' % (
        'method', 'steps', 'age[Gyr]', 's/step', 'module_s', 'retries', 'modsteps', 'E_err', 'Teff', 'Y_env', 'Y_cen', 'termination'))
    for r in results:
        if r.get('models', 0) == 0:
            print('%-26s %s' % (r['method'], r.get('status', r.get('termination', 'no history'))))
            continue
        print('%-26s %7d %9.4f %9.3f %9s %8s %7d %9.1e %8.2f %8.4f %8.4f  %s' % (
            r['method'], r['models_run'], r['age_end'] / 1e9, r['wall_per_model_s'],
            ('%.3f' % r['module_wall_s_mean']) if r.get('module_wall_s_mean') is not None else '-', r.get('retries'),
            r['module_steps'], r.get('rel_cum_E_err', float('nan')), r['Teff_end'],
            r['Y_env_end'] if r['Y_env_end'] is not None else float('nan'),
            r['Y_center_end'] if r['Y_center_end'] is not None else float('nan'), r['termination']))
    # differences at equal age relative to the first method
    ref = results[0] if results and results[0].get('models', 0) > 0 else None
    if ref is not None and len(results) > 1:
        ages = ref['history']['age']
        samples = [a for a in ages[::max(1, len(ages) // 8)]]
        print('\nTeff and Y_env relative to %s at equal age:' % ref['method'])
        print('%10s' % 'age[Gyr]' + ''.join('%22s' % r['method'][:21] for r in results[1:]))
        for a in samples:
            line = '%10.4f' % (a / 1e9)
            for r in results[1:]:
                if r.get('models', 0) == 0 or a > r['history']['age'][-1]:
                    line += '%22s' % '-'
                    continue
                dT = interp(r['history']['age'], r['history']['Teff'], a) - interp(ages, ref['history']['Teff'], a)
                dY = (interp(r['history']['age'], r['history']['Y_env'], a) - interp(ages, ref['history']['Y_env'], a)) \
                    if r['history']['Y_env'] and ref['history']['Y_env'] else float('nan')
                line += '   dTeff %6.2f dY %7.4f' % (dT, dY)
            print(line)
    for r in results:
        r.pop('history', None)
    with open(os.path.join(opts['out'], 'results.json'), 'w') as f:
        json.dump(results, f, indent=1)
    failed = [r['method'] for r in results if r.get('models', 0) == 0 or
              (opts['max_age'] is not None and r.get('age_end', 0) < 0.999 * opts['max_age'] and r.get('termination') != 'max_model_number')]
    if failed:
        print('\nruns that did not reach the requested age:', ', '.join(failed))
        if not opts['allow_failure']:
            sys.exit(1)


if __name__ == '__main__':
    main()
