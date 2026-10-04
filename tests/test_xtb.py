"""Tests of the frequency check used in the xTB step (no xTB executable required)."""
import json

from nifrec import nifrec_xtb


def _write_json(path, freqs, energy=-5.07):
    path.write_text(json.dumps({'vibrational frequencies / rcm': freqs, 'total energy': energy}))
    return str(path)


def test_small_imaginary_frequencies_below_threshold_are_ignored(tmp_path):
    path = _write_json(tmp_path / 'a.json', [0.0] * 6 + [-3.0, 150.0, 1600.0])
    assert nifrec_xtb.check_vibrational_frequencies_and_energy(path, 5.0) == (True, -5.07)


def test_significant_imaginary_frequency_is_detected(tmp_path):
    path = _write_json(tmp_path / 'b.json', [0.0] * 6 + [-10.0, 150.0, 1600.0])
    assert nifrec_xtb.check_vibrational_frequencies_and_energy(path, 5.0) == (False, -5.07)
    assert nifrec_xtb.check_vibrational_frequencies_and_energy(path, 20.0) == (True, -5.07)


def test_missing_or_broken_json(tmp_path):
    assert nifrec_xtb.check_vibrational_frequencies_and_energy(None) == (False, None)
    broken = tmp_path / 'broken.json'
    broken.write_text('{not json')
    assert nifrec_xtb.check_vibrational_frequencies_and_energy(str(broken)) == (False, None)
