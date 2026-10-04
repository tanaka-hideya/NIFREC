"""Helpers shared by the NIFREC tests.

Gaussian cannot be run on CI (commercial license). The tests of
nifrec_gaussian_optfreq therefore replace two functions with a simulated
backend (FakeGaussian):
  - nifrec_gaussian_optfreq.run_gaussian: the call of the Gaussian executable, and
  - cclib.io.ccread: the parsing of the resulting log file.
Everything else (stage logic, displacement vectors, input files, file
management, and the output .csv file) is exercised as in production.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from nifrec import nifrec_gaussian_optfreq as gof

ATOMIC_NUMBERS = {'H': 1, 'C': 6, 'N': 7, 'O': 8}

TEST_SYMBOLS = ['O', 'H', 'H']
TEST_COORDS = np.array([[0.0, 0.0, 0.1173],
                        [0.0, 0.7572, -0.4692],
                        [0.0, -0.7572, -0.4692]])
# Cartesian displacement vectors of the three normal modes of the three-atom test molecule.
# The largest-magnitude component of modes 0 and 1 is negative, so the sign convention flips them.
MODE_DISPS = np.array([
    [[0.0, 0.0, -0.60], [0.0, 0.30, 0.30], [0.0, -0.30, 0.30]],
    [[0.0, 0.0, 0.07], [0.0, -0.42, -0.56], [0.0, 0.42, -0.56]],
    [[0.0, 0.07, 0.0], [0.0, 0.56, 0.42], [0.0, -0.56, 0.42]],
])
TEXT_COLUMNS = ['imag_freqs_s0_per_cm', 'imag_freqs_s1_per_cm', 'imag_freqs_s2_per_cm']


def read_gjf_coords(path):
    """Element symbols and coordinates of a .gjf file written by NIFREC (link0/route, title, charge/coordinates)."""
    blocks = Path(path).read_text().strip().split('\n\n')
    lines = blocks[2].splitlines()[1:]
    return [line.split()[0] for line in lines], np.array([line.split()[1:4] for line in lines], dtype=float)


class FakeGaussian:
    """Simulated Gaussian backend.

    scenarios maps a job name (gname) to a list of outcomes that are consumed in call order;
    a job name without a scenario always succeeds without imaginary frequencies. Outcome keys:
      ok (default True): False simulates an abnormal termination (non-zero exit status).
      freqs (default [100, 200, 300]): vibrational frequencies in cm^-1 (negative = imaginary).
      energy (default -2000.0): final energy in eV (cclib units).
      parse_error (default False): True simulates a log file that cannot be parsed.
      move (default 0): array added to the starting geometry to mimic the optimization.
    The starting geometry is read from the .gjf file (Stage 0 and Stage 2) or taken from the
    previous job (Stage 1, Geom=AllCheck, which also requires the %oldchk file to exist).
    """

    def __init__(self, scenarios=None):
        self.scenarios = {k: list(v) for k, v in (scenarios or {}).items()}
        self.calls = []
        self._last = {}

    def run_gaussian(self, gname, gcmd):
        gjf_lines = Path(f'{gname}.gjf').read_text().splitlines()
        route = next(line for line in gjf_lines if line.startswith('#'))
        outcome = self.scenarios[gname].pop(0) if gname in self.scenarios else {}
        if 'Geom=AllCheck' in route:
            oldchk = next(line for line in gjf_lines if line.startswith('%oldchk='))[len('%oldchk='):]
            assert Path(oldchk).exists(), f'%oldchk file {oldchk} does not exist'
            symbols, coords = self._last[gname]
        else:
            symbols, coords = read_gjf_coords(f'{gname}.gjf')
        self.calls.append({'gname': gname, 'route': route, 'coords_in': np.array(coords)})
        coords_out = np.asarray(coords, dtype=float) + np.asarray(outcome.get('move', 0.0))
        Path(f'{gname}.chk').write_text('')
        Path(f'{gname}.log').write_text(json.dumps({'symbols': list(symbols),
                                                    'coords': coords_out.tolist(),
                                                    'freqs': outcome.get('freqs', [100.0, 200.0, 300.0]),
                                                    'energy': outcome.get('energy', -2000.0),
                                                    'parse_error': outcome.get('parse_error', False)}))
        if not outcome.get('ok', True):
            return False
        self._last[gname] = (list(symbols), coords_out)
        return True

    def ccread(self, logpath):
        d = json.loads(Path(logpath).read_text())
        if d['parse_error']:
            raise ValueError('simulated parse error')
        natoms = len(d['symbols'])
        nfreq = len(d['freqs'])
        disps = MODE_DISPS[:nfreq] if natoms == 3 else np.random.default_rng(natoms).normal(size=(nfreq, natoms, 3))
        return SimpleNamespace(vibfreqs=np.asarray(d['freqs'], dtype=float),
                               vibdisps=np.asarray(disps, dtype=float),
                               atomcoords=np.array([d['coords']], dtype=float),
                               atomnos=np.array([ATOMIC_NUMBERS[s] for s in d['symbols']]),
                               scfenergies=np.array([d['energy']], dtype=float))

    def stage2_inputs(self, gname):
        """Starting coordinates written to the .gjf files of Stage 2 (after Stage 0 and, if run, Stage 1)."""
        calls = [c for c in self.calls if c['gname'] == gname]
        n_skip = 1 + int(any('RCFC' in c['route'] for c in calls))
        return [c['coords_in'] for c in calls[n_skip:]]


def write_input_gjf(path, symbols=TEST_SYMBOLS, coords=TEST_COORDS, charge=0, multiplicity=1, route='#p HF/STO-3G opt'):
    lines = ['%nprocshared=1', '%chk=original.chk', route, '', 'test molecule', '', f'{charge} {multiplicity}']
    lines += [f'{s} {x:.6f} {y:.6f} {z:.6f}' for s, (x, y, z) in zip(symbols, coords)]
    Path(path).write_text('\n'.join(lines) + '\n\n')


def read_stats(path):
    return pd.read_csv(path, index_col=0, dtype={col: str for col in TEXT_COLUMNS})


def run_gjf_mode(tmp_path, names, suffix='test', infile_gaussian_recalc=None, **kwargs):
    """Run nifrec_gaussian_optfreq on .gjf inputs (one test molecule per name); return the output folder and stats."""
    infd = tmp_path / 'inputs'
    infd.mkdir(exist_ok=True)
    for name in names:
        write_input_gjf(infd / f'{name}.gjf')
    outfd = tmp_path / f'out_{suffix}'
    outfd.mkdir()
    gof.process_rows_for_goptfreq(str(outfd), None, 'xtbopt_emin_xyz', 'xTB_stats_Emin.csv', infile_gaussian_recalc, suffix,
                                  'HF/STO-3G', None, None, '=noraman', 1, 1, infd_gjf=str(infd), **kwargs)
    return outfd, read_stats(outfd / f'gaussian_{suffix}_stats.csv')


def expected_vector(mode_indices, reverse=False):
    """Expected displacement vector built from the sign-fixed MODE_DISPS."""
    vec = np.sum([gof.canonicalize_sign(MODE_DISPS[i]) for i in mode_indices], axis=0)
    vec = vec / np.linalg.norm(vec)
    return -vec if reverse else vec
