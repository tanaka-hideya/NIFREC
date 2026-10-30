"""
Tests of the Gaussian step (nifrec-gaussian-optfreq) and of the parser (nifrec-gaussian-parse).

The tests run the installed console scripts on the Gaussian input files in tests/data with PM6
(--nproc 1 --mem 3) and therefore require Gaussian 16; they take about 10 minutes in total.
The expected values below were obtained with Gaussian 16, Revision C.01.

Input files (tests/data). The file name without extension is the identifier of the molecule in the CSV files.
    0_1.gjf       NC=O: no imaginary frequency at Stage 0
    0_1_c1m2.gjf  0_1 with charge 1 and multiplicity 2: an imaginary frequency at Stage 0, removed at Stage 1
    0_1_c0m2.gjf  0_1 with charge 0 and multiplicity 2 (impossible combination): abnormal termination at Stage 0
    8_1.gjf       C[N+](=O)[O-]: imaginary frequencies remain after Stage 2 (success_disploop = -1)
    26_1.gjf      C1CCOC1: an imaginary frequency at Stage 0, removed at Stage 1
    33_1.gjf      CC(=O)c1ccccc1: an imaginary frequency at Stage 0 and at Stage 1, removed at Stage 2 (success_disploop = 1)
    284_1.gjf     CCC(=O)OC: an imaginary frequency at Stage 0, abnormal termination at Stage 1
    364_1.gjf     Cc1cc(C)cc(O)c1: two imaginary frequencies at Stage 0, removed at Stage 1
"""

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from nifrec import __version__

DATA = Path(__file__).parent / 'data'
STATS = 'gaussian_optfreq_stats.csv'
PARSE = 'gaussian_optfreq_parse.csv'

# Options of nifrec-gaussian-optfreq for each run (the recalculation reads the CSV file of 'normal')
RUNS = {'normal': [],
        'recalc': ['--option-opt-fc', 'CalcFC', '--infile-gaussian-recalc', f'normal/{STATS}'],
        'skip_stage1': ['--skip-stage1'],
        'skip_stage1_reverse_disp': ['--skip-stage1', '--reverse-disp']}

STAGE_COLUMNS = ['success_stage', 'success_disploop', 'fail_stage', 'n_imag_s0', 'n_imag_s1', 'n_imag_s2']
# Expected values of STAGE_COLUMNS for each run and molecule (None: empty cell)
EXPECTED = {
    'normal': {
        '0_1':      (0, None, None, 0, None, None),
        '0_1_c1m2': (1, None, None, 1, 0, None),
        '0_1_c0m2': (None, None, 0, None, None, None),
        '8_1':      (None, -1, None, 1, 1, 1),
        '26_1':     (1, None, None, 1, 0, None),
        '33_1':     (2, 1, None, 1, 1, 0),
        '284_1':    (None, None, 1, 1, None, None),
        '364_1':    (1, None, None, 2, 0, None),
    },
    'recalc': {  # the molecules with confid = 0 in 'normal', recalculated with opt=(CalcFC)
        '0_1_c0m2': (None, None, 0, None, None, None),
        '8_1':      (0, None, None, 0, None, None),
        '284_1':    (0, None, None, 0, None, None),
    },
    'skip_stage1': {
        '0_1':      (0, None, None, 0, None, None),
        '0_1_c1m2': (2, 0, None, 1, None, 0),
        '0_1_c0m2': (None, None, 0, None, None, None),
        '8_1':      (None, -1, None, 1, None, 1),
        '26_1':     (2, 0, None, 1, None, 0),
        '33_1':     (2, 2, None, 1, None, 0),
        '284_1':    (2, 4, None, 1, None, 0),
        '364_1':    (None, -1, None, 2, None, 2),
    },
    'skip_stage1_reverse_disp': {
        '0_1':      (0, None, None, 0, None, None),
        '0_1_c1m2': (2, 0, None, 1, None, 0),
        '0_1_c0m2': (None, None, 0, None, None, None),
        '8_1':      (None, -1, None, 1, None, 1),
        '26_1':     (2, 0, None, 1, None, 0),
        '33_1':     (2, 1, None, 1, None, 0),
        '284_1':    (2, 4, None, 1, None, 0),
        '364_1':    (2, 2, None, 2, None, 0),  # resolved only with the reversed displacement
    },
}
CHARGE_MULT = {'0_1_c1m2': (1, 2), '0_1_c0m2': (0, 2)}  # the other molecules: (0, 1)

COLUMNS = ['smiles', 'molid', 'confid', 'charge', 'multiplicity', 'total_energy_xTB', 'filepath', 'success_stage', 'success_disploop',
           'fail_stage', 'n_imag_s0', 'imag_freqs_s0_per_cm', 'energy_cclib_attr_s0', 'n_imag_s1', 'imag_freqs_s1_per_cm', 'energy_cclib_attr_s1',
           'n_imag_s2', 'imag_freqs_s2_per_cm', 'energy_cclib_attr_s2', 'rmsd_in_s0_angstrom', 'rmsd_s0_s1_angstrom', 'dE_s0_s1_kJ_per_mol',
           'rmsd_s0_s2_angstrom', 'dE_s0_s2_kJ_per_mol', 'wall_time_seconds']
ENERGY_COLUMNS = ['Final_Electronic_Energy_hartree', 'Zero_point_correction_hartree', 'Sum_of_electronic_and_zero_point_Energies_hartree',
                  'Sum_of_electronic_and_thermal_Enthalpies_hartree', 'Sum_of_electronic_and_thermal_Free_Energies_hartree']
PARSE_COLUMNS = ['functional', 'basis_set', *ENERGY_COLUMNS, 'Temperature_K', 'HOMO_eV', 'LUMO_eV']


@pytest.fixture(scope='module')
def outdir(tmp_path_factory):
    """Run nifrec-gaussian-optfreq for every entry of RUNS and nifrec-gaussian-parse on 'normal'; return the parent folder."""
    assert shutil.which('g16'), 'Gaussian 16 (g16) is required for the tests'
    root = tmp_path_factory.mktemp('gaussian')
    for run, options in RUNS.items():
        subprocess.run(['nifrec-gaussian-optfreq', '--outfolder-gaussian', run, '--infolder-gjf', str(DATA),
                        '--theory-level', 'PM6', '--nproc', '1', '--mem', '3', *options], cwd=root, check=True)
    subprocess.run(['nifrec-gaussian-parse', '--infolder-gaussian', 'normal', '--infolder-gaussian-log', 'normal/gaussian_log_optfreq',
                    '--infile', STATS, '--outfile', PARSE, '--no-homo-lumo'], cwd=root, check=True)
    return root


def read(path):
    # Output CSV file; the imaginary frequencies are kept as text ('f1;f2;...')
    return pd.read_csv(path, index_col=0, dtype={f'imag_freqs_s{stage}_per_cm': str for stage in range(3)})


def ints(row, columns):
    return tuple(None if pd.isna(row[column]) else int(row[column]) for column in columns)


def text(row, column):
    return None if pd.isna(row[column]) else row[column]


def files(folder):
    return sorted(path.name for path in folder.iterdir())


@pytest.mark.parametrize('run', RUNS)
def test_stats(outdir, run):
    df = read(outdir / run / STATS)
    assert list(df.columns) == COLUMNS
    assert sorted(df.index) == sorted(EXPECTED[run])
    for name, expected in EXPECTED[run].items():
        row = df.loc[name]
        success, n_imag = expected[0] is not None, expected[3:]
        assert ints(row, STAGE_COLUMNS) == expected
        assert ints(row, ['confid', 'charge', 'multiplicity']) == (int(success), *CHARGE_MULT.get(name, (0, 1)))
        assert row['molid'] == name and text(row, 'filepath') == (f'g_{name}.log' if success else None)
        assert pd.isna(row['smiles']) and pd.isna(row['total_energy_xTB'])  # not available for Gaussian input files
        for stage in range(3):  # the columns of a stage are filled only if the stage was run (n_imag is not None)
            filled = n_imag[stage] is not None
            freqs = text(row, f'imag_freqs_s{stage}_per_cm')
            assert text(row, f'energy_cclib_attr_s{stage}') == ('scfenergies' if filled else None)
            assert (len(freqs.split(';')) if freqs else 0) == (n_imag[stage] or 0)
            columns = ['rmsd_in_s0_angstrom'] if stage == 0 else [f'rmsd_s0_s{stage}_angstrom', f'dE_s0_s{stage}_kJ_per_mol']
            assert list(row[columns].notna()) == [filled] * len(columns)
            if filled:
                assert row[columns[0]] >= 0  # RMSD
        assert row['wall_time_seconds'] > 0


def test_output_files(outdir):
    out = outdir / 'normal'
    success = [name for name, expected in EXPECTED['normal'].items() if expected[0] is not None]
    for folder, ext in (('gjf', '.gjf'), ('log', '.log'), ('chk', '.chk')):
        assert files(out / f'gaussian_{folder}_optfreq') == sorted(f'g_{name}{ext}' for name in success)
    # imagf: Stage 0 and Stage 1 files of the molecules with imaginary frequencies and every Stage 2 trial with imaginary frequencies
    imagf = {'0_1_c1m2': ['0'], '26_1': ['0'], '284_1': ['0'], '364_1': ['0'],
             '33_1': ['0', '1', '2_0'], '8_1': ['0', '1', '2_0', '2_1', '2_2', '2_3', '2_4']}
    assert files(out / 'gaussian_imagf_optfreq') == sorted(f'g_{name}_{suffix}{ext}' for name, suffixes in imagf.items()
                                                           for suffix in suffixes for ext in ('.gjf', '.log', '.chk'))
    # working: only the files of the jobs that terminated abnormally remain
    working = files(out / 'gaussian_working_optfreq')
    assert [name for name in working if name.endswith('.gjf')] == ['g_0_1_c0m2.gjf', 'g_284_1.gjf']
    assert {name.split('.')[0] for name in working} == {'g_0_1_c0m2', 'g_284_1'}


def test_gaussian_inputs(outdir):
    # Stage 0: Link 0 and route from the options, title from the file name, charge, multiplicity, and coordinates from the input file
    lines = (outdir / 'normal' / 'gaussian_imagf_optfreq' / 'g_0_1_c1m2_0.gjf').read_text().splitlines()
    assert lines[:8] == ['%nprocshared=1', '%mem=3GB', '%chk=g_0_1_c1m2.chk', '#p PM6 opt freq=noraman', '', 'input: 0_1_c1m2.gjf', '', '1 2']
    assert [line.split()[0] for line in lines[8:14]] == ['N', 'C', 'O', 'H', 'H', 'H']
    # Stage 1: restart from the checkpoint file of Stage 0
    lines = (outdir / 'normal' / 'gaussian_gjf_optfreq' / 'g_26_1.gjf').read_text().splitlines()
    assert lines[2].startswith('%oldchk=') and lines[2].endswith('g_26_1_0.chk')
    assert lines[3:5] == ['%chk=g_26_1.chk', '#p PM6 opt=RCFC freq=noraman Guess=Read Geom=AllCheck']
    # Recalculation with --option-opt-fc CalcFC
    lines = (outdir / 'recalc' / 'gaussian_gjf_optfreq' / 'g_8_1.gjf').read_text().splitlines()
    assert lines[3] == '#p PM6 opt=(CalcFC) freq=noraman'


def test_recalc(outdir):
    normal = read(outdir / 'normal' / STATS)
    recalc = read(outdir / 'recalc' / STATS)
    assert sorted(recalc.index) == sorted(normal.index[normal['confid'] == 0])


def test_parse(outdir):
    stats = read(outdir / 'normal' / STATS)
    parse = read(outdir / 'normal' / PARSE)
    assert list(parse.columns) == COLUMNS + PARSE_COLUMNS
    assert sorted(parse.index) == sorted(stats.index[stats['confid'] != 0])
    assert (parse['confid'] == 1).all()  # every log file was parsed successfully
    pd.testing.assert_frame_equal(parse[COLUMNS], stats.loc[parse.index], check_dtype=False)  # all input columns are carried over
    assert (parse['functional'] == 'PM6').all() and (parse['basis_set'] == 'VSTO-6G').all() and (parse['Temperature_K'] == 298.15).all()
    assert parse[ENERGY_COLUMNS].notna().all().all()
    assert parse[['HOMO_eV', 'LUMO_eV']].isna().all().all()  # --no-homo-lumo


def test_log(outdir):
    log = set((outdir / 'skip_stage1_reverse_disp' / 'log_gaussian.txt').read_text().splitlines())
    assert {'theory-level: PM6', 'nproc: 1', 'mem: 3', 'disp-vec: largest', 'reverse-disp: True', 'skip-stage1: True',
            f'nifrec-version: {__version__}', 'Finish'} <= log
