"""
Example showing proper use of underscore conventions in fluoropy.
"""

import fluoropy
import numpy as np

class ImprovedPlate(fluoropy.Plate):
    """Extended plate class showing proper underscore usage."""
    
    def __init__(self, plate_format="96", name=None):
        super().__init__(plate_format, name)
        self._validation_enabled = True  # Internal setting
        self.__internal_cache = {}       # Very private cache
    
    # PUBLIC API - No underscores, in __all__
    def add_well(self, well):
        """Public method for adding wells."""
        if self._validation_enabled:
            self._validate_well(well)  # Call internal helper
        super().add_well(well)
    
    def get_statistics(self):
        """Public method for getting plate statistics."""
        return self._calculate_stats()  # Call internal helper
    
    # INTERNAL HELPERS - Single underscore
    def _validate_well(self, well):
        """Internal helper to validate well data."""
        if not well.position:
            raise ValueError("Well must have a position")
        if well.fluorescence is not None and well.fluorescence < 0:
            raise ValueError("Fluorescence cannot be negative")
    
    def _calculate_stats(self):
        """Internal helper to calculate statistics."""
        if not self.wells:
            return {"n_wells": 0}
        
        fluorescence_values = [
            w.fluorescence for w in self.wells.values() 
            if w.fluorescence is not None
        ]
        
        return {
            "n_wells": len(self.wells),
            "n_with_data": len(fluorescence_values),
            "mean_fluorescence": np.mean(fluorescence_values) if fluorescence_values else 0,
            "std_fluorescence": np.std(fluorescence_values) if fluorescence_values else 0
        }
    
    # VERY PRIVATE - Double underscore (name mangling)
    def __update_cache(self, key, value):
        """Very private method - Python mangles the name."""
        self.__internal_cache[key] = value
    
    def clear_cache(self):
        """Public method that uses private method."""
        self.__update_cache("cleared", True)
        self.__internal_cache.clear()


# Demonstrate the conventions
def show_underscore_conventions():
    """Show how underscore conventions work in practice."""
    
    plate = ImprovedPlate("96", "Demo Plate")
    
    # PUBLIC API - Users should use these
    well = fluoropy.Well("A1", fluorescence=1000)
    plate.add_well(well)                    # ✅ Public method
    stats = plate.get_statistics()          # ✅ Public method
    print("Stats:", stats)
    
    # INTERNAL API - Users can access but shouldn't
    plate._validate_well(well)              # ⚠️ Works but shouldn't use
    internal_stats = plate._calculate_stats()  # ⚠️ Works but internal
    
    # VERY PRIVATE - Name mangling makes it harder to access
    try:
        plate.__update_cache("test", "value")  # ❌ AttributeError!
    except AttributeError as e:
        print(f"Can't access __private method: {e}")
    
    # But you can still access if you know the mangled name
    plate._ImprovedPlate__update_cache("test", "value")  # ⚠️ Possible but very bad
    
    print("Cache contents:", plate._ImprovedPlate__internal_cache)


# Show what should be in __all__ for this class
__all__ = [
    "ImprovedPlate",           # ✅ Public class
    "show_underscore_conventions"  # ✅ Public function
    # NOT included: _validate_well, _calculate_stats, __update_cache
]


if __name__ == "__main__":
    show_underscore_conventions()
