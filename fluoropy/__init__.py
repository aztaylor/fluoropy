"""
fluoropy: A Python package for fluorescence assay calculations and utilities to
aid in experimental design and analysis.
"""

__version__ = "0.1.0"
__author__ = "Aleczander Taylor"
__email__ = "aztaylor76@fastmail.com"

# Import analysis modules (users can access as fluoropy.analysis.*)
# Import utility modules (users can access as fluoropy.utils.*)
from . import analysis, utils

# Import main modules/classes here for easy access
from .core import (  # noqa: F401  (plotting is exposed as fluoropy.plotting)
    Fluorophore,
    Plate,
    Sample,
    SampleFrame,
    Well,
    plotting,
)
from .core.fluorophore import fluorophore_db

# Make key classes/functions available at package level
__all__ = [
    "Plate",
    "Well",
    "Sample",
    "SampleFrame",
    "Fluorophore",
    "fluorophore_db",
    "utils",
    "analysis",
]
