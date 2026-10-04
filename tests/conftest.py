import pytest

from nifrec import nifrec_gaussian_optfreq as gof
from nifrec_testutils import FakeGaussian


@pytest.fixture
def fake_gaussian(monkeypatch):
    """Install a simulated Gaussian backend; call it with a scenario dict (see FakeGaussian)."""
    def install(scenarios=None):
        fake = FakeGaussian(scenarios)
        monkeypatch.setattr(gof, 'run_gaussian', fake.run_gaussian)
        monkeypatch.setattr(gof.cclib.io, 'ccread', fake.ccread)
        return fake
    return install
