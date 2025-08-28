"""
Example showing good vs bad API design for fluoropy.
"""

# =============================================================================
# GOOD API DESIGN - What users should be able to do easily
# =============================================================================

import fluoropy

def good_api_examples():
    """Examples of a well-designed, user-friendly API."""
    
    # 1. SIMPLE IMPORTS - Everything they need from one place
    plate = fluoropy.Plate("96")  # Not: from fluoropy.core.plate import Plate
    well = fluoropy.Well("A1", 1000)
    assay = fluoropy.EndpointAssay("test", "FITC")
    
    # 2. SENSIBLE DEFAULTS - Works with minimal configuration
    plate = fluoropy.Plate()  # Defaults to 96-well
    well = fluoropy.Well("A1", 1000)  # Minimal required info
    
    # 3. CONSISTENT NAMING - Similar operations have similar names
    plate.add_well(well)        # add_*
    plate.get_well("A1")        # get_*
    plate.remove_well("A1")     # remove_* (if you add this)
    
    # 4. FLUENT INTERFACE - Methods can be chained
    experiment = (fluoropy.PlateSet("Experiment")
                  .add_plate(plate1)
                  .add_plate(plate2)
                  .combine_plates())
    
    # 5. HELPFUL UTILITIES - Common tasks are easy
    combined = fluoropy.combine_plates([plate1, plate2])  # One function call
    
    # 6. CLEAR RETURN TYPES - Users know what they get back
    results = assay.analyze()  # Always returns Dict
    df = plate.get_fluorescence_data()  # Always returns DataFrame


# =============================================================================
# BAD API DESIGN - What to avoid
# =============================================================================

def bad_api_examples():
    """Examples of confusing, hard-to-use API design."""
    
    # 1. COMPLEX IMPORTS - Users need to know internal structure
    from fluoropy.core.plate import Plate
    from fluoropy.core.wells import Well  
    from fluoropy.analysis.endpoint import EndpointAssay
    
    # 2. TOO MANY REQUIRED PARAMETERS
    plate = Plate(
        format="96", 
        validate=True, 
        allow_duplicates=False,
        auto_sort=True,
        metadata_validation="strict"
    )
    
    # 3. INCONSISTENT NAMING
    plate.addWell(well)         # camelCase
    plate.get_well_data("A1")   # snake_case
    plate.DeleteWell("A1")      # PascalCase
    
    # 4. UNCLEAR RETURN TYPES
    result = assay.analyze()    # Sometimes Dict, sometimes List, sometimes None?
    
    # 5. NO HELPFUL SHORTCUTS
    # Users have to write complex code for simple tasks
    experiment = PlateSet()
    for plate in plates:
        experiment.add_plate(plate)
    combined = experiment.combine_plates("average")
    # Instead of: combined = combine_plates(plates, method="average")


# =============================================================================
# YOUR CURRENT API EVALUATION
# =============================================================================

def evaluate_fluoropy_api():
    """Your current API is actually quite good! Here's why:"""
    
    # ✅ GOOD: Simple imports
    import fluoropy
    plate = fluoropy.Plate("96")
    
    # ✅ GOOD: Sensible defaults
    plate = fluoropy.Plate()  # Defaults to 96-well
    
    # ✅ GOOD: Consistent naming
    plate.add_well(well)
    plate.get_well("A1")
    
    # ✅ GOOD: Clear purpose classes
    endpoint_assay = fluoropy.EndpointAssay(name, fluorophore, plate)
    
    # ✅ GOOD: Utility functions
    combined = fluoropy.combine_plates(plates)
    
    # 🚀 SUGGESTIONS FOR IMPROVEMENT:
    
    # 1. Add more convenience constructors
    def suggested_improvements():
        # Current: Verbose well creation
        well = fluoropy.Well("A1", fluorescence=1000, concentration=10, compound="X")
        
        # Suggested: Convenience methods
        plate = fluoropy.Plate.from_excel("data.xlsx")
        plate = fluoropy.Plate.from_csv("data.csv")
        well = fluoropy.Well.sample("A1", 1000, 10, "X")  # clear purpose
        well = fluoropy.Well.control("H1", 100)
        well = fluoropy.Well.blank("H12", 50)
        
        # 2. Method chaining
        results = (fluoropy.EndpointAssay("test", "FITC")
                  .add_plate(plate)
                  .subtract_background()
                  .analyze())


if __name__ == "__main__":
    print("API Design Examples - see source code for details!")
    good_api_examples()
    evaluate_fluoropy_api()
