"""
Unit and regression test for the lj package.
"""

# Import package, test suite, and other packages as needed
import lj
import pytest
import sys

def test_lj_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "lj" in sys.modules
