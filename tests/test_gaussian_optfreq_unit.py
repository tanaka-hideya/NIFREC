"""Unit tests of the building blocks of nifrec_gaussian_optfreq (no Gaussian required)."""
import io
import itertools
import os
import sys
import textwrap
from types import SimpleNamespace

import numpy as np
import pytest

from nifrec import nifrec_gaussian_optfreq as gof
from nifrec_testutils import MODE_DISPS, TEST_COORDS


# ---------- sign convention and displacement vector ----------

def test_canonicalize_sign_makes_largest_component_positive():
    vec = np.array([[0.1, -0.5, 0.2], [0.3, 0.0, -0.1]])
    assert gof.canonicalize_sign(vec) == pytest.approx(-vec)
    assert gof.canonicalize_sign(-vec) == pytest.approx(-vec)


def test_canonicalize_sign_uses_first_component_for_ties():
    assert gof.canonicalize_sign(np.array([0.4, -0.4])) == pytest.approx([0.4, -0.4])
    assert gof.canonicalize_sign(np.array([-0.4, 0.4])) == pytest.approx([0.4, -0.4])


def test_displacement_vector_is_independent_of_the_signs_of_the_input_modes():
    freqs = np.array([-50.0, -20.0, 300.0])
    for disp_vec in ('sum', 'largest'):
        reference = gof.combined_imag_vector(freqs, MODE_DISPS, disp_vec, False)
        for signs in itertools.product([1.0, -1.0], repeat=3):
            disps = MODE_DISPS * np.array(signs)[:, None, None]
            assert gof.combined_imag_vector(freqs, disps, disp_vec, False) == pytest.approx(reference)


def test_displacement_vector_sum_and_largest():
    freqs = np.array([-20.0, -50.0, 300.0])
    expected_sum = gof.canonicalize_sign(MODE_DISPS[0]) + gof.canonicalize_sign(MODE_DISPS[1])
    assert (gof.combined_imag_vector(freqs, MODE_DISPS, 'sum', False)
            == pytest.approx(expected_sum / np.linalg.norm(expected_sum)))
    expected_largest = gof.canonicalize_sign(MODE_DISPS[1])  # most negative frequency: mode 1
    assert (gof.combined_imag_vector(freqs, MODE_DISPS, 'largest', False)
            == pytest.approx(expected_largest / np.linalg.norm(expected_largest)))


def test_displacement_vector_single_mode_and_reverse():
    freqs = np.array([-30.0, 200.0, 300.0])
    vec_sum = gof.combined_imag_vector(freqs, MODE_DISPS, 'sum', False)
    assert vec_sum == pytest.approx(gof.combined_imag_vector(freqs, MODE_DISPS, 'largest', False))
    assert np.linalg.norm(vec_sum) == pytest.approx(1.0)
    assert gof.combined_imag_vector(freqs, MODE_DISPS, 'sum', True) == pytest.approx(-vec_sum)
    assert gof.combined_imag_vector(freqs, MODE_DISPS, 'largest', True) == pytest.approx(-vec_sum)


def test_summarize_imag():
    assert gof.summarize_imag(np.array([-20.48123, 5.0, -100.0, 300.0])) == (2, '-100.0000;-20.4812')
    assert gof.summarize_imag(np.array([1.0, 2.0])) == (0, '')


# ---------- energies, RMSD, and log parsing ----------

def test_final_energy_uses_the_highest_level_available():
    assert gof.final_energy(SimpleNamespace(scfenergies=np.array([-2000.5, -2001.0]))) == pytest.approx(-2001.0)
    data = SimpleNamespace(scfenergies=np.array([-2001.0]), mpenergies=np.array([[-2002.0, -2003.0]]))
    assert gof.final_energy(data) == pytest.approx(-2003.0)
    data.ccenergies = np.array([-2004.0])
    assert gof.final_energy(data) == pytest.approx(-2004.0)


def test_parse_freq_and_disp(monkeypatch):
    data = SimpleNamespace(vibfreqs=np.array([-30.0, 200.0, 300.0]), vibdisps=MODE_DISPS,
                           atomcoords=np.array([TEST_COORDS + 0.1, TEST_COORDS]), scfenergies=np.array([-2000.5, -2001.0]))
    monkeypatch.setattr(gof.cclib.io, 'ccread', lambda path: data)
    freqs, disps, coords, energy = gof.parse_freq_and_disp('dummy.log')
    assert freqs == pytest.approx([-30.0, 200.0, 300.0])
    assert disps == pytest.approx(MODE_DISPS)
    assert coords == pytest.approx(TEST_COORDS)  # last geometry
    assert energy == pytest.approx(-2001.0)


def test_rmsd_aligned_known_value():
    a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    b = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    assert gof.rmsd_aligned(a, b) == pytest.approx(0.5)


def test_rmsd_aligned_is_invariant_to_rotation_and_translation():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(8, 3))
    rot, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(rot) < 0:
        rot[:, 0] *= -1
    moved = coords @ rot.T + np.array([1.0, -2.0, 3.0])
    assert gof.rmsd_aligned(coords, moved) == pytest.approx(0.0, abs=1e-8)


def test_rmsd_aligned_does_not_superpose_mirror_images():
    chiral = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    assert gof.rmsd_aligned(chiral, chiral * np.array([1.0, 1.0, -1.0])) > 0.1


# ---------- input files and Gaussian execution ----------

def test_write_gjf(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    gof.write_gjf(2, 4, 'g_mol', '#p HF/STO-3G opt freq=noraman', 'smiles: O', 0, 1, 'O 0.0 0.0 0.0\n')
    assert (tmp_path / 'g_mol.gjf').read_text() == (
        '%nprocshared=2\n%mem=4GB\n%chk=g_mol.chk\n#p HF/STO-3G opt freq=noraman\n\nsmiles: O\n\n0 1\nO 0.0 0.0 0.0\n\n')
    gof.write_gjf(2, 4, 'g_rcfc', '#p HF/STO-3G opt=RCFC freq Guess=Read Geom=AllCheck', '', '', '', '', '/x/g_0.chk')
    assert (tmp_path / 'g_rcfc.gjf').read_text() == (
        '%nprocshared=2\n%mem=4GB\n%oldchk=/x/g_0.chk\n%chk=g_rcfc.chk\n'
        '#p HF/STO-3G opt=RCFC freq Guess=Read Geom=AllCheck\n\n\n')


def _write_script(path, body):
    path.write_text('#!/bin/sh\n' + body)
    os.chmod(path, 0o755)
    return str(path)


@pytest.mark.skipif(sys.platform.startswith('win'), reason='POSIX shell scripts are used as fake Gaussian executables')
def test_run_gaussian_with_fake_executables(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'g_mol.gjf').write_text('dummy\n')
    writes_log = _write_script(tmp_path / 'g_log.sh', 'echo done > "${1%.gjf}.log"\n')
    writes_out = _write_script(tmp_path / 'g_out.sh', 'echo done > "${1%.gjf}.out"\n')
    fails = _write_script(tmp_path / 'g_fail.sh', 'echo error > "${1%.gjf}.log"\nexit 1\n')
    silent = _write_script(tmp_path / 'g_silent.sh', 'exit 0\n')

    assert gof.run_gaussian('g_mol', writes_log)
    (tmp_path / 'g_mol.log').unlink()
    assert gof.run_gaussian('g_mol', writes_out)
    assert (tmp_path / 'g_mol.log').exists() and not (tmp_path / 'g_mol.out').exists()
    (tmp_path / 'g_mol.log').unlink()
    assert not gof.run_gaussian('g_mol', fails)
    (tmp_path / 'g_mol.log').unlink()
    assert not gof.run_gaussian('g_mol', silent)


# ---------- command-line interface ----------

def test_cli_defaults_and_validation():
    args = gof._parse_cli_args(['--outfolder-gaussian', 'out', '--infolder-gjf', 'inputs'])
    assert args.infolder_gjf == 'inputs' and args.infolder_xtb is None
    assert (args.disp_vec, args.reverse_disp, args.skip_stage1) == ('sum', False, False)
    assert (args.base_disp, args.max_repeat) == (0.1, 5)
    args = gof._parse_cli_args(['--outfolder-gaussian', 'out', '--infolder-xtb', 'xtb', '--disp-vec', 'largest',
                                '--reverse-disp', '--skip-stage1'])
    assert (args.disp_vec, args.reverse_disp, args.skip_stage1) == ('largest', True, True)
    invalid = [['--outfolder-gaussian', 'out'],
               ['--outfolder-gaussian', 'out', '--infolder-xtb', 'x', '--infolder-gjf', 'y'],
               ['--outfolder-gaussian', 'out', '--infolder-gjf', 'y', '--disp-vec', 'mean'],
               ['--outfolder-gaussian', 'out', '--infolder-gjf', 'y', '--imag-vec']]
    for argv in invalid:
        with pytest.raises(SystemExit):
            gof._parse_cli_args(argv)


def test_help_shows_the_module_description(monkeypatch):
    out = io.StringIO()
    monkeypatch.setattr(sys, 'stdout', out)
    with pytest.raises(SystemExit):
        gof._parse_cli_args(['--help'])
    help_text = out.getvalue()
    description = textwrap.dedent(gof.__doc__.split('Description:', 1)[1]).strip()
    for line in description.splitlines():  # the description of --help is identical to the docstring
        assert line.rstrip() in help_text
    for col in gof.INFO_COLUMNS:  # every added output column is documented
        assert col in help_text
