"""
fluoropy: A Python package for fluorescence assay calculations and utilities to
aid in experimental design and analysis.
"""

__version__ = "0.1.0"
__author__ = "Aleczander Taylor"
__email__ = "aztaylor76@fastmail.com"

# Import main modules/classes here for easy access
from .core import Plate, Well, Sample, SampleFrame, Fluorophore, Assay
from .core.fluorophore import fluorophore_db
from .core.assay import EndpointAssay

# Import utility modules (users can access as fluoropy.utils.*)
from . import utils

# Import analysis modules (users can access as fluoropy.analysis.*)
from . import analysis

# Make key classes/functions available at package level
__all__ = [
    "Plate",
    "Well",
    "Sample",
    "SampleFrame",
    "Fluorophore",
    "fluorophore_db",
    "Assay",
    "EndpointAssay",
    "utils",
    "analysis",
]
