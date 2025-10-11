#!/usr/bin/env python3

# Simple test without imports to debug
print("Starting simple test...")

class SimpleWell:
    def __init__(self, well_id):
        self.well_id = well_id
        self.concentration = None
        self._excluded = False
        self.measurements = {}

    def exclude(self):
        self._excluded = True

    def is_excluded(self):
        return self._excluded

    def add_measurement(self, name, data):
        self.measurements[name] = data

class SimpleSample:
    def __init__(self):
        self.wells = []

    def add_well(self, well):
        self.wells.append(well)

    def get_concentrations(self):
        concentrations = []
        seen = set()

        for well in self.wells:
            if (well.concentration is not None and
                not well.is_excluded() and
                well.concentration not in seen):
                concentrations.append(well.concentration)
                seen.add(well.concentration)

        return concentrations

# Test the concentration ordering
wells = []
for i, conc in enumerate([0.1, 0.2, 0.3, 0.4]):
    well = SimpleWell(f'A{i+1}')
    well.concentration = conc
    wells.append(well)

sample = SimpleSample()
for well in wells:
    sample.add_well(well)

print("Before exclusion:", sample.get_concentrations())

# Exclude first well
wells[0].exclude()
print("After excluding first well:", sample.get_concentrations())

# The expected result should be [0.2, 0.3, 0.4] NOT [0.2, 0.3, 0.4, 0.1]
expected = [0.2, 0.3, 0.4]
actual = sample.get_concentrations()

if actual == expected:
    print("✅ Concentration ordering is correct")
else:
    print(f"❌ Issue found: expected {expected}, got {actual}")

print("Test completed.")