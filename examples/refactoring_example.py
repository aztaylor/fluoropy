"""
Example: How to refactor analysis methods from Plate class

This shows the BEFORE and AFTER of moving analysis logic
to separate modules while keeping the Plate class clean.
"""

# ======================================================================
# BEFORE: Analysis mixed with data structure
# ======================================================================

class PlateOld:
    def __init__(self):
        self.wells = {}

    def calculate_cv_old_way(self, well_list):
        """This was INSIDE the Plate class - mixing concerns!"""
        values = [self.wells[pos].fluorescence for pos in well_list]
        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)
        return (std_val / mean_val) * 100


# ======================================================================
# AFTER: Clean separation of concerns
# ======================================================================

# 1. Plate class focuses ONLY on data structure
class Plate:
    def __init__(self):
        self.wells = {}

    def get_well(self, position):
        """Data access - this belongs in Plate"""
        return self.wells.get(position)

    def get_fluorescence_data(self):
        """Data access - this belongs in Plate"""
        # Keep this method - it's data access, not analysis
        pass


# 2. Analysis functions in separate module
from fluoropy.analysis.statistics import calculate_cv

# Usage becomes:
plate = Plate()
cv_result = calculate_cv(plate, ['A1', 'A2', 'A3'])

# ======================================================================
# 🎯 KEY BENEFITS of this approach:
# ======================================================================

"""
✅ Single Responsibility:
   - Plate = data storage + basic access
   - Analysis modules = computations

✅ Testability:
   - Test data structures separately
   - Test analysis functions separately

✅ Maintainability:
   - Easy to add new analysis without changing Plate
   - Analysis functions can be reused across different data types

✅ API Clarity:
   - Users import specific analysis functions they need
   - Plate class stays focused and simple
"""
