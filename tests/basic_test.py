#!/usr/bin/env python3

import sys
import os

# Change to the fluoropy directory
os.chdir('/Users/alec/Documents/SideProjects/fluoropy')
sys.path.insert(0, '.')

print("=== Basic Import Test ===")

# Test each import step by step
try:
    print("1. Importing numpy...")
    import numpy as np
    print("   ✅ numpy imported successfully")

    print("2. Importing pandas...")
    import pandas as pd
    print("   ✅ pandas imported successfully")

    print("3. Importing fluoropy package...")
    import fluoropy
    print("   ✅ fluoropy package imported")

    print("4. Importing fluoropy.core module...")
    import fluoropy.core
    print("   ✅ fluoropy.core module imported")

    print("5. Importing Well class...")
    from fluoropy.core.well import Well
    print("   ✅ Well class imported")

    print("6. Importing Plate class...")
    from fluoropy.core.plate_new import Plate
    print("   ✅ Plate class imported")

    print("7. Importing Sample class...")
    from fluoropy.core.sample import Sample
    print("   ✅ Sample class imported")

    print("8. Importing SampleFrame class...")
    from fluoropy.core.sampleframe_new import SampleFrame
    print("   ✅ SampleFrame class imported")

    print("\n=== Creating Basic Objects ===")

    print("9. Creating a Well...")
    well = Well("A1", 0, 0)
    print(f"   ✅ Well created: {well}")

    print("10. Creating a Plate...")
    plate = Plate("96", "test_plate")
    print(f"    ✅ Plate created: {plate}")

    print("11. Accessing well from plate...")
    well_from_plate = plate['A1']
    print(f"    ✅ Well A1 from plate: {well_from_plate}")

    print("\n🎉 ALL BASIC TESTS PASSED!")
    print("The new architecture imports and basic functionality work correctly.")

except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
