"""
Test layout calculation for large number of samples
"""
import sys
sys.path.insert(0, '/Users/alec/Documents/SideProjects/fluoropy')

import numpy as np
from fluoropy.core.plate import Plate
from fluoropy.core.sampleframe import SampleFrame


def test_layout_calculation():
    """Test that layout calculation works for large sample counts"""
    print("="*70)
    print("Testing Layout Calculation for Various Sample Counts")
    print("="*70)

    test_cases = [
        (2, 2),     # 4 subplots
        (3, 4),     # 12 subplots
        (5, 4),     # 20 subplots
        (10, 4),    # 40 subplots
        (15, 4),    # 60 subplots
        (20, 4),    # 80 subplots
    ]

    for n_samples, n_conc in test_cases:
        print(f"\n{n_samples} samples × {n_conc} concentrations = {n_samples * n_conc} subplots")

        # Create plate with adequate size (384-well has 16 rows × 24 cols = 384 wells)
        plate = Plate(plate_format="384", name="test")

        # Fill with samples
        row, col = 0, 0
        well_count = 0
        for s in range(n_samples):
            sample_name = f's{s}'
            for c in range(n_conc):
                # Stop if we run out of wells
                if row >= plate.rows or col >= plate.cols:
                    break

                conc = 0.1 * (c + 1)
                well = plate.get_well_by_position(row, col)

                if well is None:
                    print(f"  ✗ Could not get well at ({row}, {col})")
                    break

                well.set_sample_info(sample_name, conc)
                well.time_points = np.linspace(0, 24, 15)
                well.time_series['OD600'] = np.random.random(15)

                col += 1
                well_count += 1
                if col >= plate.cols:
                    col = 0
                    row += 1

        # Create SampleFrame and plot
        frame = SampleFrame(plate)

        try:
            fig, axes = frame.plot_replicate_time_series(
                'OD600',
                show_mean=False,  # Faster without mean
            )

            n_subplots = len(axes)
            figsize = fig.get_size_inches()

            # Calculate grid from figsize
            expected_cols = int(round(figsize[0] / 3.5))
            expected_rows = int(round(figsize[1] / 3.5))

            aspect_ratio = figsize[0] / figsize[1]

            print(f"  ✓ Figure size: {figsize[0]:.1f}\" × {figsize[1]:.1f}\"")
            print(f"  ✓ Grid: ~{expected_cols} cols × {expected_rows} rows")
            print(f"  ✓ Aspect ratio: {aspect_ratio:.2f} (golden ratio ≈ 1.4)")
            print(f"  ✓ Subplots created: {n_subplots}")

            if aspect_ratio < 1.2 or aspect_ratio > 1.6:
                print(f"  ⚠ Aspect ratio slightly outside golden ratio")

            # Clean up
            import matplotlib.pyplot as plt
            plt.close(fig)

        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            return False

    print("\n" + "="*70)
    print("✓ All layout tests PASSED!")
    print("="*70)
    return True


if __name__ == "__main__":
    success = test_layout_calculation()
    exit(0 if success else 1)
