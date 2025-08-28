"""Tests for fluoropy package."""
import pytest
from fluoropy import __version__


def test_version():
    """Test that version is defined."""
    assert __version__ == "0.1.0"


def test_import():
    """Test that the package can be imported."""
    import fluoropy
    assert fluoropy is not None
