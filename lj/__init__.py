"""
LJ
The goal is to reproduce Rahman's paper (Correlations in the Motion of Atoms in Liquid Argon), with 864 particles interacting with a Lennard-Jones potential and compute corresponding properties
"""

# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Add imports here
from .gr import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
