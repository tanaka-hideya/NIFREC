"""NIFREC: A Gaussian-Integrated Automated Workflow for Local-Minimum Geometry Optimization with No Imaginary Frequencies"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('nifrec')
except PackageNotFoundError:  # e.g., running from a source tree without installation
    __version__ = 'unknown'
