"""Run the repository's numerical reproduction, with optional figure builds."""
from pathlib import Path
import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures', action='store_true')
    parser.add_argument('--supplement', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    env = os.environ.copy()
    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[name] = '1'
    env['PYTHONIOENCODING'] = 'utf-8'
    env['MPLBACKEND'] = 'Agg'
    scripts = ['recompute_statistics.py', 'verify_statistics.py']
    if args.figures:
        scripts.append('build_figures.py')
    if args.supplement:
        scripts.extend(['build_supplement.py', 'prepare_workbook_data.py'])
    for script in scripts:
        print('Running ' + script, flush=True)
        subprocess.run([sys.executable, str(root / '05_REPRODUCIBILITY' / script)],
                       cwd=root, env=env, check=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
