#!/usr/bin/env python3
"""
Test script to demonstrate the concentration ordering issue when excluding wells.
"""

from fluoropy.core.well import Well
from fluoropy.core.sample import Sample

def test_concentration_ordering():
    # Create wells with different concentrations
    wells = []
    concentrations = [0.1, 0.2, 0.3, 0.4]

    for i, conc in enumerate(concentrations):
        well = Well(f'A{i+1}')
        well.concentration = conc
        well.add_time_series('measurement', [1.0 * (i+1), 2.0 * (i+1), 3.0 * (i+1)])
        wells.append(well)

    print('Original wells and concentrations:')
    for well in wells:
        print(f'  {well.well_id}: concentration={well.concentration}, data={well.get_measurement("measurement")}')

    # Create sample and add wells
    sample = Sample('test')
    for well in wells:
        sample.add_well(well)

    print('\nBefore exclusion:')
    print('Concentrations from get_concentrations():', sample.get_concentrations())

    sample.calculate_statistics(['measurement'])
    print('Time series shape:', sample.time_series['measurement'].shape)
    print('Time series data:\n', sample.time_series['measurement'])
    print('Sample.concentrations:', sample.concentrations)

    # Now exclude the first well (concentration 0.1)
    print(f'\n--- Excluding first well ({wells[0].well_id}, conc={wells[0].concentration}) ---')
    wells[0].exclude()

    print('\nAfter exclusion:')
    print('Concentrations from get_concentrations():', sample.get_concentrations())

    sample.calculate_statistics(['measurement'])
    print('Time series shape:', sample.time_series['measurement'].shape)
    print('Time series data:\n', sample.time_series['measurement'])
    print('Sample.concentrations:', sample.concentrations)

    # The issue: concentration 0.1 should be completely gone, but it might be appearing at the end

if __name__ == "__main__":
    test_concentration_ordering()