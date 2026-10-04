"""Tests of nifrec_gaussian_parse with simulated cclib data (no Gaussian required)."""
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from nifrec import nifrec_gaussian_optfreq as gof
from nifrec import nifrec_gaussian_parse as gparse


def _cclib_data(freqs):
    return SimpleNamespace(metadata={'success': True, 'functional': 'HF', 'basis_set': 'STO-3G'},
                           optdone=[1], vibfreqs=np.array(freqs), scfenergies=np.array([-2000.0, -2001.0]),
                           zpve=0.02, enthalpy=-73.5, freeenergy=-73.53, temperature=298.15,
                           homos=np.array([4]), moenergies=[np.array([-20.0, -15.0, -10.0, -8.0, -5.0, 2.0, 4.0])])


def test_parse_carries_over_all_columns(tmp_path, monkeypatch):
    logfd = tmp_path / 'logs'
    logfd.mkdir()
    for name in ('g_a.log', 'g_b.log', 'g_d.log'):
        (logfd / name).write_text('dummy')
    info = dict.fromkeys(gof.INFO_COLUMNS)
    base = {'smiles': 'O', 'molid': 'a', 'confid': 1, 'charge': 0, 'multiplicity': 1, 'total_energy_xTB': None,
            'filepath': None, 'success_stage': None, 'success_disploop': None}
    status = {
        'a': {**base, **info, 'filepath': 'g_a.log', 'success_stage': 0, 'n_imag_s0': 0, 'imag_freqs_s0_per_cm': '',
              'rmsd_in_s0_angstrom': 0.01, 'wall_time_seconds': 1.5},
        'b': {**base, **info, 'molid': 'b', 'filepath': 'g_b.log', 'success_stage': 2, 'success_disploop': 1,
              'n_imag_s0': 2, 'imag_freqs_s0_per_cm': '-50.0000;-20.0000', 'n_imag_s1': 1, 'imag_freqs_s1_per_cm': '-40.0000',
              'n_imag_s2': 0, 'imag_freqs_s2_per_cm': '', 'rmsd_s0_s2_angstrom': 0.05, 'dE_s0_s2_kJ_per_mol': -0.25},
        'c': {**base, **info, 'molid': 'c', 'confid': 0, 'fail_stage': 1, 'n_imag_s0': 1, 'imag_freqs_s0_per_cm': '-30.0000'},
        'd': {**base, **info, 'molid': 'd', 'filepath': 'g_d.log', 'success_stage': 0, 'n_imag_s0': 0},
    }
    stats = pd.DataFrame.from_dict(status, orient='index')
    for col in ('success_stage', 'success_disploop', 'fail_stage', 'n_imag_s0', 'n_imag_s1', 'n_imag_s2'):
        stats[col] = stats[col].astype('Int64')
    stats.to_csv(tmp_path / 'gaussian_test_stats.csv')

    data = {'g_a.log': _cclib_data([100.0, 200.0, 300.0]),
            'g_b.log': _cclib_data([100.0, 200.0, 300.0]),
            'g_d.log': _cclib_data([-10.0, 200.0, 300.0])}  # an imaginary frequency is rejected by the parser
    monkeypatch.setattr(gparse.cclib.io, 'ccread', lambda path: data[os.path.basename(path)])
    gparse.process_rows_for_gparse(str(tmp_path), str(logfd), 'gaussian_test_stats.csv', 'parsed.csv', False)

    raw = pd.read_csv(tmp_path / 'parsed.csv', index_col=0, dtype=str, keep_default_na=False)
    assert list(raw.index) == ['a', 'b', 'd']  # confid = 0 in the input is not parsed
    assert list(raw.columns[:9 + len(gof.INFO_COLUMNS)]) == list(stats.columns)  # all input columns are carried over
    assert raw.loc['b', 'success_stage'] == '2' and raw.loc['b', 'success_disploop'] == '1'  # kept as integers
    assert raw.loc['b', 'imag_freqs_s0_per_cm'] == '-50.0000;-20.0000' and raw.loc['b', 'imag_freqs_s1_per_cm'] == '-40.0000'
    assert raw.loc['b', 'dE_s0_s2_kJ_per_mol'] == '-0.25'
    assert raw.loc['a', 'confid'] == '1' and raw.loc['d', 'confid'] == '0'
    parsed = pd.read_csv(tmp_path / 'parsed.csv', index_col=0)
    # energies are converted back to hartree with the same factor that cclib used (27.21138505 eV/hartree)
    assert parsed.loc['a', 'Final_SCF_Energy_hartree'] == pytest.approx(-2001.0 / 27.21138505, rel=1e-12)
    assert parsed.loc['a', 'HOMO_eV'] == pytest.approx(-5.0)
    assert parsed.loc['a', 'LUMO_eV'] == pytest.approx(2.0)
    assert parsed.loc['a', 'Sum_of_electronic_and_thermal_Free_Energies_hartree'] == pytest.approx(-73.53)


def test_parse_rejects_identical_file_names(tmp_path):
    with pytest.raises(ValueError, match='identical'):
        gparse.process_rows_for_gparse(str(tmp_path), str(tmp_path), 'same.csv', 'same.csv', False)
