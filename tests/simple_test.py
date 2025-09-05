#!/usr/bin/env python3
import sys
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

print("Testing imports step by step...")

try:
    import fluoropy
    print("✅ fluoropy module imported")
except Exception as e:
    print(f"❌ fluoropy import failed: {e}")
    exit(1)

try:
    import fluoropy.core
    print("✅ fluoropy.core imported")
except Exception as e:
    print(f"❌ fluoropy.core import failed: {e}")
    exit(1)

try:
    from fluoropy.core import Plate
    print("✅ Plate imported")
except Exception as e:
    print(f"❌ Plate import failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

try:
    from fluoropy.core import Well
    print("✅ Well imported")
except Exception as e:
    print(f"❌ Well import failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

try:
    from fluoropy.core import Sample
    print("✅ Sample imported")
except Exception as e:
    print(f"❌ Sample import failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

try:
    from fluoropy.core import SampleFrame
    print("✅ SampleFrame imported")
except Exception as e:
    print(f"❌ SampleFrame import failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

print("🎉 All imports successful!")

# Quick test
plate = Plate(plate_format="96", name="test")
print(f"Created plate: {plate}")
print(f"Well A1: {plate['A1']}")
