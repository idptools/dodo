"""
Unit and regression test for the dodo package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import dodo


def test_dodo_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "dodo" in sys.modules
