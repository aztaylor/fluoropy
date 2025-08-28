"""
Module showing how to hide implementation details effectively.
"""

# This module demonstrates module-level privacy

# ============================================================================
# PRIVATE MODULE VARIABLES AND FUNCTIONS
# ============================================================================

# Private constants (not in __all__)
_DEFAULT_BACKGROUND = 50.0
_MAX_FLUORESCENCE = 100000.0
_SUPPORTED_FORMATS = ["96", "384", "1536"]

# Private helper functions
def _validate_fluorescence_value(value):
    """Internal validation function."""
    if value < 0:
        raise ValueError("Fluorescence cannot be negative")
    if value > _MAX_FLUORESCENCE:
        raise ValueError(f"Fluorescence too high (max: {_MAX_FLUORESCENCE})")
    return True

def _calculate_z_score(value, mean, std):
    """Internal statistical calculation."""
    if std == 0:
        return 0
    return (value - mean) / std

def _format_well_position(row, col):
    """Internal function to format well positions."""
    return f"{chr(ord('A') + row)}{col + 1}"

# Private class for internal use
class _InternalCache:
    """Internal caching mechanism - not for public use."""
    
    def __init__(self):
        self._cache = {}
    
    def get(self, key):
        return self._cache.get(key)
    
    def set(self, key, value):
        self._cache[key] = value
    
    def clear(self):
        self._cache.clear()

# Global private cache instance
_global_cache = _InternalCache()

# ============================================================================
# PUBLIC API
# ============================================================================

class WellValidator:
    """Public class for validating well data."""
    
    def __init__(self):
        self._cache = _InternalCache()  # Uses private class
    
    def validate_well_data(self, position, fluorescence):
        """
        Public method to validate well data.
        
        Parameters
        ----------
        position : str
            Well position (e.g., 'A1')
        fluorescence : float
            Fluorescence value
            
        Returns
        -------
        bool
            True if valid, raises exception if not
        """
        # Use private validation functions
        _validate_fluorescence_value(fluorescence)
        
        # Cache the validation result
        cache_key = f"{position}:{fluorescence}"
        self._cache.set(cache_key, True)
        
        return True
    
    def get_validation_stats(self):
        """Get statistics about validations performed."""
        # This method could use private statistical functions
        return {
            "validations_cached": len(self._cache._cache),
            "max_allowed_fluorescence": _MAX_FLUORESCENCE,
            "default_background": _DEFAULT_BACKGROUND
        }


def create_well_layout(plate_format):
    """
    Public function to create well layout.
    
    Parameters
    ----------
    plate_format : str
        Plate format ("96", "384", or "1536")
        
    Returns
    -------
    list
        List of well positions
    """
    if plate_format not in _SUPPORTED_FORMATS:  # Uses private constant
        raise ValueError(f"Unsupported format. Use: {_SUPPORTED_FORMATS}")
    
    # Use private helper function
    if plate_format == "96":
        rows, cols = 8, 12
    elif plate_format == "384":
        rows, cols = 16, 24
    else:  # 1536
        rows, cols = 32, 48
    
    positions = []
    for row in range(rows):
        for col in range(cols):
            positions.append(_format_well_position(row, col))  # Private function
    
    return positions


# ============================================================================
# DEFINE PUBLIC API
# ============================================================================

# Only these items can be imported with "from module import *"
__all__ = [
    "WellValidator",
    "create_well_layout"
]

# Notice what's NOT in __all__:
# - _validate_fluorescence_value
# - _calculate_z_score  
# - _format_well_position
# - _InternalCache
# - _global_cache
# - _DEFAULT_BACKGROUND
# - _MAX_FLUORESCENCE
# - _SUPPORTED_FORMATS


def demonstrate_module_privacy():
    """Show how module-level privacy works."""
    
    print("=== PUBLIC API USAGE ===")
    
    # Use public classes and functions
    validator = WellValidator()
    is_valid = validator.validate_well_data("A1", 1000.0)
    print(f"Validation result: {is_valid}")
    
    stats = validator.get_validation_stats()
    print(f"Validation stats: {stats}")
    
    layout = create_well_layout("96")
    print(f"96-well layout: {len(layout)} positions")
    print(f"First few positions: {layout[:5]}")
    
    print("\n=== WHAT'S HIDDEN ===")
    
    # These work from within the module but shouldn't be used externally:
    print(f"Private constant accessible here: {_DEFAULT_BACKGROUND}")
    print(f"Private function accessible here: {_format_well_position(0, 0)}")
    
    # Show what's available for import
    print(f"\nPublic API (__all__): {__all__}")
    
    # Show all module contents (including private)
    all_names = [name for name in globals() if not name.startswith('__')]
    public_names = [name for name in all_names if not name.startswith('_')]
    private_names = [name for name in all_names if name.startswith('_')]
    
    print(f"Public names: {public_names}")
    print(f"Private names: {private_names}")


if __name__ == "__main__":
    demonstrate_module_privacy()
