"""End-to-end test: RDKit conformers -> xTB -> Gaussian step (simulated Gaussian backend).

Requires RDKit, MORFEUS, and the xtb executable (all provided by environment.yml);
skipped otherwise.
"""
import shutil

import pandas as pd
import pytest

from nifrec import nifrec_gaussian_optfreq as gof
from nifrec_testutils import read_stats


def test_rdkit_xtb_gaussian_workflow(tmp_path, fake_gaussian, monkeypatch):
    pytest.importorskip('morfeus')
    if shutil.which('xtb') is None:
        pytest.skip('xtb executable not found')
    from nifrec import nifrec_rdkit, nifrec_xtb

    monkeypatch.chdir(tmp_path)
    infile = tmp_path / 'mols.csv'
    infile.write_text('name,smi\nethanol,CCO\npropane,CCC\n')

    rdkit_out = tmp_path / 'rdkit'
    rdkit_out.mkdir()
    nifrec_rdkit.process_rows_for_rdkit(str(rdkit_out), str(infile), 'rdkit_stats.csv', 3, 0.5,
                                        njobs=1, smicol='smi')
    assert pd.read_csv(rdkit_out / 'rdkit_stats.csv', index_col=0)['success_rdkit_confgen'].all()

    xtb_out = tmp_path / 'xtb'
    xtb_out.mkdir()
    nifrec_xtb.process_rows_for_xtb(str(xtb_out), str(rdkit_out), 'rdkit_stats.csv', 'xTB_stats_all.csv',
                                    'xTB_stats_Emin.csv', max_nconf=2, max_repeat=5, njobs=1)
    assert (pd.read_csv(xtb_out / 'xTB_stats_Emin.csv', index_col=0)['confid'] > 0).all()

    fake_gaussian()  # every Gaussian job succeeds at Stage 0
    gauss_out = tmp_path / 'gaussian'
    gauss_out.mkdir()
    gof.process_rows_for_goptfreq(str(gauss_out), str(xtb_out), 'xtbopt_emin_xyz', 'xTB_stats_Emin.csv', None,
                                  'it', 'HF/STO-3G', None, None, '=noraman', 1, 1)
    stats = read_stats(gauss_out / 'gaussian_it_stats.csv')
    assert list(stats.index) == ['ethanol', 'propane']
    assert (stats['success_stage'] == 0).all()
    assert (stats['charge'] == 0).all() and (stats['multiplicity'] == 1).all()
    assert (stats['rmsd_in_s0_angstrom'] < 1e-6).all()
