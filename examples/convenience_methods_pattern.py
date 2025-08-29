"""
Example: Adding convenience methods to Plate class that delegate to analysis functions.

This shows how to provide easy access while maintaining clean separation.
"""

from typing import List, Dict, Optional
import pandas as pd
import numpy as np

# Import analysis functions (the real implementation)
from ..analysis.statistics import calculate_cv as _calculate_cv
from ..analysis.statistics import calculate_z_factor as _calculate_z_factor
from ..analysis.normalization import normalize_to_controls as _normalize_to_controls
from ..analysis.normalization import percent_inhibition as _percent_inhibition

class PlateWithConvenienceMethods:
    """
    Example showing how to add convenience methods to Plate class
    that delegate to analysis functions.

    ✅ Benefits:
    - Easy to use: plate.calculate_cv(['A1', 'A2'])
    - Still maintains separation of concerns
    - Analysis functions remain standalone and testable
    - Users can choose either API style
    """

    def __init__(self, plate_format: str = "96", name: Optional[str] = None):
        self.format = plate_format
        self.name = name
        self.wells = {}
        # ... rest of init ...

    # ======================================================================
    # CONVENIENCE METHODS - Delegate to analysis functions
    # ======================================================================

    def calculate_cv(self, well_list: List[str], timepoint: Optional[int] = None) -> float:
        """
        Calculate coefficient of variation for wells.

        Convenience method that delegates to analysis.statistics.calculate_cv

        Args:
            well_list: List of well positions
            timepoint: Specific timepoint for kinetic data

        Returns:
            CV as percentage

        Example:
            >>> cv = plate.calculate_cv(['A1', 'A2', 'A3'])
        """
        return _calculate_cv(self, well_list, timepoint)

    def calculate_z_factor(self, positive_controls: List[str],
                          negative_controls: List[str],
                          timepoint: Optional[int] = None) -> float:
        """
        Calculate Z-factor for assay quality.

        Convenience method that delegates to analysis.statistics.calculate_z_factor

        Args:
            positive_controls: Positive control well positions
            negative_controls: Negative control well positions
            timepoint: Specific timepoint for kinetic data

        Returns:
            Z-factor value

        Example:
            >>> z_factor = plate.calculate_z_factor(['A1', 'A2'], ['B1', 'B2'])
        """
        return _calculate_z_factor(self, positive_controls, negative_controls, timepoint)

    def normalize_to_controls(self, test_wells: List[str],
                             positive_controls: List[str],
                             negative_controls: List[str],
                             timepoint: Optional[int] = None) -> Dict[str, float]:
        """
        Normalize test wells to control wells.

        Convenience method that delegates to analysis.normalization.normalize_to_controls

        Args:
            test_wells: Wells to normalize
            positive_controls: Positive control wells
            negative_controls: Negative control wells
            timepoint: Specific timepoint for kinetic data

        Returns:
            Dict mapping well positions to normalized values

        Example:
            >>> normalized = plate.normalize_to_controls(['C1', 'C2'], ['A1', 'A2'], ['B1', 'B2'])
        """
        return _normalize_to_controls(self, test_wells, positive_controls, negative_controls, timepoint)

    def percent_inhibition(self, test_wells: List[str],
                          control_wells: List[str],
                          timepoint: Optional[int] = None) -> Dict[str, float]:
        """
        Calculate percent inhibition relative to controls.

        Convenience method that delegates to analysis.normalization.percent_inhibition

        Args:
            test_wells: Wells with test compounds
            control_wells: Control wells (no inhibition)
            timepoint: Specific timepoint for kinetic data

        Returns:
            Dict mapping well positions to percent inhibition values

        Example:
            >>> inhibition = plate.percent_inhibition(['C1', 'C2'], ['A1', 'A2'])
        """
        return _percent_inhibition(self, test_wells, control_wells, timepoint)


# ======================================================================
# USAGE COMPARISON - Both APIs available
# ======================================================================

def usage_examples():
    """Show how both API styles work."""

    plate = PlateWithConvenienceMethods()

    # Method 1: Convenience methods (easy)
    cv = plate.calculate_cv(['A1', 'A2', 'A3'])
    z_factor = plate.calculate_z_factor(['A1', 'A2'], ['B1', 'B2'])
    normalized = plate.normalize_to_controls(['C1', 'C2'], ['A1', 'A2'], ['B1', 'B2'])

    # Method 2: Direct analysis functions (explicit)
    from fluoropy.analysis.statistics import calculate_cv, calculate_z_factor
    from fluoropy.analysis.normalization import normalize_to_controls

    cv = calculate_cv(plate, ['A1', 'A2', 'A3'])
    z_factor = calculate_z_factor(plate, ['A1', 'A2'], ['B1', 'B2'])
    normalized = normalize_to_controls(plate, ['C1', 'C2'], ['A1', 'A2'], ['B1', 'B2'])


# ======================================================================
# 🎯 DESIGN PRINCIPLES for convenience methods:
# ======================================================================

"""
✅ DO:
- Keep convenience methods thin (just delegate)
- Import analysis functions with underscore prefix (_calculate_cv)
- Maintain identical signatures to analysis functions
- Document that they're convenience methods
- Keep analysis functions as the "source of truth"

❌ DON'T:
- Duplicate analysis logic in convenience methods
- Make convenience methods do more than delegate
- Change function signatures between convenience and analysis versions
- Remove the standalone analysis functions
"""
