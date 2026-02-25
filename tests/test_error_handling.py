"""
Test error handling in plot_replicate_time_series
"""
import sys
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame


def test_error_handling():
    """Test that error handling works properly"""
    print("="*60)
    print("Testing Error Handling")
    print("="*60)

    # Test 1: Invalid sample ID
    print("\n1. Testing invalid sample ID...")
    plate = Plate(name="test")
    well = plate.get_well_by_position(0, 0)
    well.set_sample_info('s14', 0.1)
    well.time_points = None
    well.time_series['OD600'] = []

    frame = SampleFrame(plate)

    try:
        frame.plot_replicate_time_series('OD600', sample_ids=['nonexistent'])
        print("   ✗ FAILED - Should have raised ValueError")
        return False
    except ValueError as e:
        print(f"   ✓ PASSED - Caught: {str(e)[:60]}...")

    # Test 2: Invalid measurement
    print("\n2. Testing invalid measurement...")
    try:
        frame.plot_replicate_time_series('INVALID_MEASUREMENT', sample_ids=['s14'])
        print("   ✗ FAILED - Should have raised ValueError")
        return False
    except ValueError as e:
        print(f"   ✓ PASSED - Caught: {str(e)[:60]}...")

    # Test 3: Missing time_points
    print("\n3. Testing missing time_points...")
    try:
        frame.plot_replicate_time_series('OD600', sample_ids=['s14'])
        print("   ✗ FAILED - Should have raised RuntimeError")
        return False
    except RuntimeError as e:
        print(f"   ✓ PASSED - Caught: {str(e)[:60]}...")

    # Test 4: Empty sample list
    print("\n4. Testing empty sample list...")
    try:
        frame.plot_replicate_time_series('OD600', sample_ids=[])
        print("   ✗ FAILED - Should have raised ValueError")
        return False
    except ValueError as e:
        print(f"   ✓ PASSED - Caught: {str(e)[:60]}...")

    print("\n" + "="*60)
    print("✓ All error handling tests PASSED!")
    print("="*60)
    return True


if __name__ == "__main__":
    success = test_error_handling()
    exit(0 if success else 1)
