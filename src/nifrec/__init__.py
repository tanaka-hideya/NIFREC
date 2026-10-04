"""NIFREC: automated local-minimum geometry optimization with no imaginary vibrational frequencies."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('nifrec')
except PackageNotFoundError:  # e.g., running from a source tree without installation
    __version__ = 'unknown'
