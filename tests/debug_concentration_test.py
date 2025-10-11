#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.getcwd())

from fluoropy.core.well import Well
from fluoropy.core.sample import Sample

def test_concentration_ordering():
    print("Testing concentration ordering issue...")

    # Create wells with different concentrations
    wells = []
    concentrations = [0.1, 0.2, 0.3, 0.4]

    for i, conc in enumerate(concentrations):
        well = Well(f'A{i+1}')
        well.concentration = conc
        well.add_time_series('measurement', [1.0 * (i+1), 2.0 * (i+1), 3.0 * (i+1)])
        wells.append(well)

    print('\nOriginal wells:')
    for well in wells:
        print(f'  {well.well_id}: conc={well.concentration}')

    # Create sample
    sample = Sample('test')
    for well in wells:
        sample.add_well(well)

    print('\nBefore exclusion:')
    concs_before = sample.get_concentrations()
    print('Concentrations from get_concentrations():', concs_before)

    sample.calculate_statistics(['measurement'])
    print('Time series shape:', sample.time_series['measurement'].shape)
    print('Sample.concentrations:', sample.concentrations)
    print('Time series data:')
    print(sample.time_series['measurement'])

    # Exclude first well (concentration 0.1)
    print(f'\n--- Excluding first well ({wells[0].well_id}, conc={wells[0].concentration}) ---')
    wells[0].exclude()

    print('\nAfter exclusion:')
    concs_after = sample.get_concentrations()
    print('Concentrations from get_concentrations():', concs_after)

    sample.calculate_statistics(['measurement'])
    print('Time series shape after exclusion:', sample.time_series['measurement'].shape)
    print('Sample.concentrations after exclusion:', sample.concentrations)
    print('Time series data after exclusion:')
    print(sample.time_series['measurement'])

    # Check if the issue is present
    if len(concs_after) != len(concs_before) - 1:
        print(f"\n❌ ISSUE: Expected {len(concs_before) - 1} concentrations after exclusion, got {len(concs_after)}")

    if 0.1 in concs_after:
        print(f"\n❌ ISSUE: Excluded concentration 0.1 still appears in concentrations array")
    else:
        print(f"\n✅ GOOD: Excluded concentration 0.1 properly removed from concentrations")

if __name__ == "__main__":
    test_concentration_ordering()