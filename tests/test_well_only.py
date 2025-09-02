#!/usr/bin/env python3
import sys
import os
print("Starting test...")

# Just test basic imports step by step
try:
    print("Importing Well...")
    from fluoropy.core.well import Well
    print("Well imported successfully")

    print("Creating a Well...")
    w = Well("A1", 0, 0)
    print(f"Well created: {w.well_id}")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
