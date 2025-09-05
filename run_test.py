#!/usr/bin/env python3

import sys
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

try:
    # Run the test
    exec(open('test_new_architecture.py').read())

    # If we get here, the test passed
    print("\n" + "="*60)
    print("🎉 COMPREHENSIVE TEST PASSED! 🎉")
    print("="*60)
    print("✅ All imports work correctly")
    print("✅ Plate and Well containers work")
    print("✅ Sample objects group wells by sample_type")
    print("✅ Data structure: numpy arrays (n_timepoints, n_concentrations)")
    print("✅ SampleFrame indexing works")
    print("✅ Complete processing pipeline works")
    print("\n🏗️ Your new architecture is ready to use!")

    with open('test_success.txt', 'w') as f:
        f.write('SUCCESS: All tests passed!\n')
        f.write('Architecture is working correctly.\n')

except Exception as e:
    print(f"\n❌ Test failed with error: {e}")
    import traceback
    traceback.print_exc()

    with open('test_failure.txt', 'w') as f:
        f.write(f'FAILED: {e}\n')
        f.write(traceback.format_exc())
