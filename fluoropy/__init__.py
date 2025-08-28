"""
fluoropy: A Python package for fluorescence assay calculations and utilities to
aid in experimental design and analysis.
"""

__version__ = "0.1.0"
__author__ = "Aleczander Taylor"
__email__ = "aztaylor76@fastmail.com"

# Import main modules/classes here for easy access
from .core import Plate, Well, Fluorophore, Assay, PlateSet, combine_plates
from .core.fluorophore import fluorophore_db
from .core.assay import EndpointAssay

# Make key classes/functions available at package level
__all__ = [
    "Plate",
    "Well",
    "Fluorophore",
    "fluorophore_db",
    "Assay",
    "EndpointAssay",
    "PlateSet",
    "combine_plates",
]
