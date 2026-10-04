"""Optional regression test with real Gaussian output files.

Place Gaussian 'opt freq' log files (for example, small molecules computed with NIFREC,
including at least one with imaginary frequencies) in tests/data/gaussian/ to check that
they are parsed correctly with the installed cclib. The test is skipped when no log file
is present (Gaussian output files are not distributed with NIFREC).
"""
from pathlib import Path

import numpy as np
import pytest

from nifrec import nifrec_gaussian_optfreq as gof

LOG_DIR = Path(__file__).parent / 'data' / 'gaussian'


def test_real_gaussian_logs():
    logs = sorted(LOG_DIR.glob('*.log')) if LOG_DIR.is_dir() else []
    if not logs:
        pytest.skip('no Gaussian log files in tests/data/gaussian')
    for log in logs:
        freqs, disps, coords, energy = gof.parse_freq_and_disp(str(log))
        assert disps.shape == (len(freqs), len(coords), 3), log.name
        assert coords.shape[1] == 3 and np.isfinite(energy), log.name
        if gof.detect_imag(freqs):
            vec = gof.combined_imag_vector(freqs, disps, 'sum', False)
            assert np.linalg.norm(vec) == pytest.approx(1.0), log.name
