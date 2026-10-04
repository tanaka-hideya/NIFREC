"""Tests of the Stage 0-2 workflow of nifrec_gaussian_optfreq with a simulated Gaussian backend."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nifrec import nifrec_gaussian_optfreq as gof
from nifrec_testutils import TEST_COORDS, expected_vector, read_stats, run_gjf_mode, write_input_gjf

IMAG = {'freqs': [-50.0, 200.0, 300.0]}
OK = {'freqs': [100.0, 200.0, 300.0]}
EV_TO_KJMOL = 96.4853364596  # cclib.parser.utils.convertor


def _files(folder):
    return sorted(p.name for p in Path(folder).iterdir())


def test_stage0_success(tmp_path, fake_gaussian):
    fake_gaussian({'g_mol1': [OK]})
    outfd, stats = run_gjf_mode(tmp_path, ['mol1'])
    row = stats.loc['mol1']
    assert row['confid'] == 1 and row['success_stage'] == 0 and row['filepath'] == 'g_mol1.log'
    assert pd.isna(row['success_disploop']) and pd.isna(row['fail_stage'])
    assert row['n_imag_s0'] == 0 and pd.isna(row['imag_freqs_s0_per_cm'])
    assert pd.isna(row['n_imag_s1']) and pd.isna(row['n_imag_s2'])
    assert row['rmsd_in_s0_angstrom'] == pytest.approx(0.0, abs=1e-6)
    assert pd.isna(row['rmsd_s0_s1_angstrom']) and pd.isna(row['dE_s0_s1_kJ_per_mol'])
    assert row['wall_time_seconds'] >= 0
    assert (row['charge'], row['multiplicity']) == (0, 1)
    assert _files(outfd / 'gaussian_gjf_test') == ['g_mol1.gjf']
    assert _files(outfd / 'gaussian_log_test') == ['g_mol1.log']
    assert _files(outfd / 'gaussian_working_test') == []
    gjf_text = (outfd / 'gaussian_gjf_test' / 'g_mol1.gjf').read_text()
    assert '#p HF/STO-3G opt freq=noraman\n\ninput: mol1.gjf\n\n0 1\n' in gjf_text


def test_stage1_success_records_rcfc_restart(tmp_path, fake_gaussian):
    move = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.05], [0.0, 0.0, 0.0]])
    fake = fake_gaussian({'g_mol1': [{'freqs': [-50.0, 200.0, 300.0], 'energy': -2000.0},
                                     {'freqs': [10.0, 200.0, 300.0], 'energy': -2000.001, 'move': move}]})
    outfd, stats = run_gjf_mode(tmp_path, ['mol1'])
    row = stats.loc['mol1']
    assert row['success_stage'] == 1 and pd.isna(row['success_disploop'])
    assert row['n_imag_s0'] == 1 and row['imag_freqs_s0_per_cm'] == '-50.0000'
    assert row['n_imag_s1'] == 0 and pd.isna(row['imag_freqs_s1_per_cm'])
    assert row['dE_s0_s1_kJ_per_mol'] == pytest.approx(-0.001 * EV_TO_KJMOL, abs=1e-6)  # rounded to 1e-6 kJ/mol
    assert row['rmsd_s0_s1_angstrom'] == pytest.approx(gof.rmsd_aligned(TEST_COORDS, TEST_COORDS + move), abs=1e-6)
    assert row['rmsd_s0_s1_angstrom'] > 0
    assert 'opt=RCFC freq=noraman Guess=Read Geom=AllCheck' in fake.calls[1]['route']
    assert _files(outfd / 'gaussian_imagf_test') == ['g_mol1_0.chk', 'g_mol1_0.gjf', 'g_mol1_0.log']
    assert _files(outfd / 'gaussian_working_test') == []


def test_stage2_displacements_start_from_the_same_structure(tmp_path, fake_gaussian):
    fake = fake_gaussian({'g_mol1': [IMAG, {'freqs': [-40.0, 200.0, 300.0]}, {'freqs': [-30.0, 200.0, 300.0]}, OK]})
    outfd, stats = run_gjf_mode(tmp_path, ['mol1'])
    row = stats.loc['mol1']
    assert (row['success_stage'], row['success_disploop']) == (2, 1) and pd.isna(row['fail_stage'])
    assert row['imag_freqs_s1_per_cm'] == '-40.0000'
    assert row['n_imag_s2'] == 0 and pd.isna(row['imag_freqs_s2_per_cm'])
    assert row['rmsd_s0_s2_angstrom'] > 0
    starts = fake.stage2_inputs('g_mol1')
    assert len(starts) == 2
    for k, coords in enumerate(starts):  # trial k: the Stage 1 structure displaced by 0.1 * (k + 1)
        assert coords == pytest.approx(TEST_COORDS + expected_vector([0]) * 0.1 * (k + 1), abs=2e-6)
    assert _files(outfd / 'gaussian_imagf_test') == sorted(
        f'g_mol1_{s}.{e}' for s in ('0', '1', '2_0') for e in ('chk', 'gjf', 'log'))
    assert _files(outfd / 'gaussian_working_test') == []


def test_stage2_vector_options_and_skip_stage1(tmp_path, fake_gaussian):
    two_imag = {'freqs': [-50.0, -20.0, 300.0]}
    settings = [('sum', False, [0, 1]), ('largest', False, [0]), ('sum', True, [0, 1]), ('largest', True, [0])]
    for n, (disp_vec, reverse, modes) in enumerate(settings):
        fake = fake_gaussian({'g_mol1': [two_imag, OK]})
        _, stats = run_gjf_mode(tmp_path, ['mol1'], suffix=f'opt{n}', base_disp=0.15,
                                disp_vec=disp_vec, reverse_disp=reverse, skip_stage1=True)
        row = stats.loc['mol1']
        assert (row['success_stage'], row['success_disploop']) == (2, 0)
        assert pd.isna(row['n_imag_s1'])  # Stage 1 was skipped
        assert row['imag_freqs_s0_per_cm'] == '-50.0000;-20.0000'
        assert not any('RCFC' in c['route'] for c in fake.calls)
        start = fake.stage2_inputs('g_mol1')[0]
        assert start == pytest.approx(TEST_COORDS + expected_vector(modes, reverse) * 0.15, abs=2e-6)


def test_unresolved_imaginary_frequencies(tmp_path, fake_gaussian):
    fake_gaussian({'g_mol1': [IMAG, IMAG, {'freqs': [-30.0, 200.0, 300.0]}, {'freqs': [-12.34, 200.0, 300.0]}]})
    outfd, stats = run_gjf_mode(tmp_path, ['mol1'], max_repeat=2)
    row = stats.loc['mol1']
    assert row['confid'] == 0 and row['success_disploop'] == -1
    assert pd.isna(row['success_stage']) and pd.isna(row['fail_stage']) and pd.isna(row['filepath'])
    assert row['n_imag_s2'] == 1 and row['imag_freqs_s2_per_cm'] == '-12.3400'  # last trial
    assert not pd.isna(row['rmsd_s0_s2_angstrom']) and not pd.isna(row['dE_s0_s2_kJ_per_mol'])
    assert 'g_mol1_2_1.log' in _files(outfd / 'gaussian_imagf_test')
    assert _files(outfd / 'gaussian_working_test') == []


def test_abnormal_termination_keeps_files_and_continues(tmp_path, fake_gaussian):
    fake_gaussian({'g_mol1': [IMAG, {'ok': False}], 'g_mol2': [OK]})
    outfd, stats = run_gjf_mode(tmp_path, ['mol1', 'mol2'])
    failed = stats.loc['mol1']
    assert failed['confid'] == 0 and failed['fail_stage'] == 1
    assert pd.isna(failed['success_stage']) and pd.isna(failed['success_disploop'])
    assert failed['n_imag_s0'] == 1 and pd.isna(failed['n_imag_s1'])
    assert _files(outfd / 'gaussian_working_test') == ['g_mol1.chk', 'g_mol1.gjf', 'g_mol1.log']
    assert stats.loc['mol2', 'success_stage'] == 0


def test_stage2_failure_clears_stage2_columns(tmp_path, fake_gaussian):
    fake_gaussian({'g_mol1': [IMAG, IMAG, {'freqs': [-30.0, 200.0, 300.0]}, {'ok': False}]})
    _, stats = run_gjf_mode(tmp_path, ['mol1'])
    row = stats.loc['mol1']
    assert row['fail_stage'] == 2 and row['confid'] == 0 and pd.isna(row['success_disploop'])
    assert pd.isna(row['n_imag_s2']) and pd.isna(row['rmsd_s0_s2_angstrom']) and pd.isna(row['dE_s0_s2_kJ_per_mol'])
    assert row['n_imag_s1'] == 1


def test_errors_are_recorded_and_the_batch_continues(tmp_path, fake_gaussian, monkeypatch):
    fake_gaussian({'g_mol1': [{'parse_error': True}], 'g_mol2': [IMAG], 'g_mol3': [OK]})

    def broken_vector(*args, **kwargs):
        raise RuntimeError('simulated bug')
    monkeypatch.setattr(gof, 'combined_imag_vector', broken_vector)
    _, stats = run_gjf_mode(tmp_path, ['mol1', 'mol2', 'mol3'], skip_stage1=True)
    assert stats.loc['mol1', 'fail_stage'] == 0 and pd.isna(stats.loc['mol1', 'n_imag_s0'])
    assert stats.loc['mol2', 'fail_stage'] == 2 and stats.loc['mol2', 'confid'] == 0
    assert stats.loc['mol3', 'success_stage'] == 0


def test_recalculation_of_failed_molecules(tmp_path, fake_gaussian):
    fake_gaussian({'g_001': [{'ok': False}], 'g_002': [OK]})
    outfd, _ = run_gjf_mode(tmp_path, ['001', '002'], suffix='first')
    fake_gaussian({'g_001': [OK]})
    outfd2, _ = run_gjf_mode(tmp_path, ['001', '002'], suffix='second',
                             infile_gaussian_recalc=str(outfd / 'gaussian_first_stats.csv'))
    stats2 = pd.read_csv(outfd2 / 'gaussian_second_stats.csv', index_col=0, dtype=str)
    assert list(stats2.index) == ['001']  # identifiers with leading zeros are matched as text
    assert stats2.loc['001', 'success_stage'] == '0'


def test_gjf_input_errors(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # the working directory is restored after the test
    cases = {
        'duplicate': ({'a.gjf': None, 'a.com': None}, 'unique'),
        'whitespace': ({'a b.gjf': None}, 'whitespace'),
        'no_charge': ({'a.gjf': '%oldchk=x.chk\n#p HF/STO-3G Geom=AllCheck Guess=Read\n\n'}, 'could not be read'),
        'bohr': ({'a.gjf': 'BOHR'}, 'angstrom'),
        'empty': ({'notes.txt': 'not an input file'}, 'No .gjf or .com'),
    }
    for name, (files, message) in cases.items():
        infd = tmp_path / name
        infd.mkdir()
        for fname, text in files.items():
            if text is None:
                write_input_gjf(infd / fname)
            elif text == 'BOHR':
                write_input_gjf(infd / fname, route='#p HF/STO-3G Units=Bohr')
            else:
                (infd / fname).write_text(text)
        outfd = tmp_path / f'out_{name}'
        outfd.mkdir()
        with pytest.raises(ValueError, match=message):
            gof.process_rows_for_goptfreq(str(outfd), None, 'xtbopt_emin_xyz', 'xTB_stats_Emin.csv', None, 'test',
                                          'HF/STO-3G', None, None, '=noraman', 1, 1, infd_gjf=str(infd))


def test_charge_and_multiplicity_from_gjf(tmp_path, fake_gaussian):
    fake_gaussian()
    infd = tmp_path / 'inputs'
    infd.mkdir()
    write_input_gjf(infd / 'cation.gjf', symbols=['N', 'H', 'H'], charge=1, multiplicity=2)
    outfd = tmp_path / 'out'
    outfd.mkdir()
    gof.process_rows_for_goptfreq(str(outfd), None, 'xtbopt_emin_xyz', 'xTB_stats_Emin.csv', None, 'test',
                                  'UHF/STO-3G', None, None, '=noraman', 1, 1, infd_gjf=str(infd))
    stats = read_stats(outfd / 'gaussian_test_stats.csv')
    assert (stats.loc['cation', 'charge'], stats.loc['cation', 'multiplicity']) == (1, 2)
    assert 'input: cation.gjf\n\n1 2\nN ' in (outfd / 'gaussian_gjf_test' / 'g_cation.gjf').read_text()


def test_main_writes_log(tmp_path, fake_gaussian, monkeypatch):
    fake_gaussian({'g_mol1': [OK]})
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'stdout', sys.stdout)  # main() redirects and closes sys.stdout; restored after the test
    (tmp_path / 'inputs').mkdir()
    write_input_gjf(tmp_path / 'inputs' / 'mol1.gjf')
    gof.main(['--outfolder-gaussian', 'out', '--infolder-gjf', 'inputs', '--theory-level', 'PM6', '--suffix', 'pm6', '--nproc', '1'])
    log = (tmp_path / 'out' / 'log_gaussian.txt').read_text()
    for text in ('nifrec-version:', 'hostname:', 'cpu:', 'disp-vec: sum', 'skip-stage1: False', 'Finish'):
        assert text in log
    stats = read_stats(tmp_path / 'out' / 'gaussian_pm6_stats.csv')
    assert list(stats.columns) == (['smiles', 'molid', 'confid', 'charge', 'multiplicity', 'total_energy_xTB', 'filepath',
                                    'success_stage', 'success_disploop'] + gof.INFO_COLUMNS)


def test_xtb_csv_input_with_charge_and_multiplicity_from_smiles(tmp_path, fake_gaussian):
    pytest.importorskip('rdkit')
    xtbfd = tmp_path / 'xtb'
    (xtbfd / 'xtbopt_emin_xyz').mkdir(parents=True)
    rows = {0: ('O', 1, ['O', 'H', 'H']), 1: ('[CH2]', 1, ['C', 'H', 'H']), 2: ('CC', 0, None)}
    records = {}
    for number, (smiles, confid, symbols) in rows.items():
        filepath = f'xTB_{number}_1.xyz' if confid else ''
        records[number] = {'smiles': smiles, 'molid': number, 'confid': confid, 'total_energy_xTB': -5.0,
                           'filepath': filepath}
        if symbols:
            lines = ['3', ' energy: -5.0'] + [f'{s} {x:.8f} {y:.8f} {z:.8f}' for s, (x, y, z) in zip(symbols, TEST_COORDS)]
            (xtbfd / 'xtbopt_emin_xyz' / filepath).write_text('\n'.join(lines) + '\n')
    pd.DataFrame.from_dict(records, orient='index').to_csv(xtbfd / 'xTB_stats_Emin.csv')
    fake_gaussian()
    outfd = tmp_path / 'gauss'
    outfd.mkdir()
    gof.process_rows_for_goptfreq(str(outfd), str(xtbfd), 'xtbopt_emin_xyz', 'xTB_stats_Emin.csv', None, 'csv',
                                  'HF/STO-3G', None, None, '=noraman', 1, 1)
    stats = read_stats(outfd / 'gaussian_csv_stats.csv')
    assert list(stats.index) == [0, 1]  # confid = 0 in the xTB summary is skipped
    assert list(stats['multiplicity']) == [1, 3]  # number of radical electrons + 1 ([CH2]: high spin)
    assert list(stats['success_stage']) == [0, 0]
    assert (stats['rmsd_in_s0_angstrom'] < 1e-6).all()
    gjf_text = (outfd / 'gaussian_gjf_csv' / 'g_1_1.gjf').read_text()
    assert 'smiles: [CH2]\n\n0 3\nC ' in gjf_text
