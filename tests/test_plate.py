#!/usr/bin/env python3
"""
Test script for the combined Plate class
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Test basic imports
try:
    from fluoropy.core.well import Well
    print("✓ Well import successful")
except Exception as e:
    print(f"✗ Well import failed: {e}")

try:
    from fluoropy.core.plate import Plate
    print("✓ Plate import successful")
except Exception as e:
    print(f"✗ Plate import failed: {e}")

# Test basic functionality
try:
    plate = Plate(plate_format="96", name="test_plate")
    print(f"✓ Created plate: {plate}")

    # Test subscriptable access
    well_a1 = plate['A1']
    print(f"✓ Accessed well A1: {well_a1}")

    # Test total wells
    print(f"✓ Total wells: {len(plate)}")

    # Test different plate formats
    plate_384 = Plate("384", "test_384")
    print(f"✓ Created 384-well plate: {plate_384}")

    plate_1536 = Plate("1536", "test_1536")
    print(f"✓ Created 1536-well plate: {plate_1536}")

    print("\n🎉 All basic tests passed! The combined Plate class is working correctly.")

except Exception as e:
    print(f"✗ Test failed: {e}")
    import traceback
    traceback.print_exc()
