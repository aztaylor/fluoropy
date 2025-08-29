"""
Complete example showing the new analysis architecture.

This demonstrates how analysis functions are now separated from
the Plate class, following the Single Responsibility Principle.
"""

import fluoropy as fp
import numpy as np

def main():
    print("🧬 Fluoropy Analysis Architecture Demo")
    print("="*50)

    # ======================================================================
    # 1. Create plate with data (Plate class focused on DATA STRUCTURE)
    # ======================================================================

    plate = fp.Plate(plate_format="96", name="Test Plate")

    # Add some sample wells
    for row in range(3):  # A, B, C
        for col in range(4):  # 1, 2, 3, 4
            well_pos = f"{chr(65 + row)}{col + 1}"

            # Simulate different well types
            if row == 0:  # Row A = positive controls
                well_type = "positive_control"
                fluorescence = np.random.normal(1000, 50)
            elif row == 1:  # Row B = negative controls
                well_type = "negative_control"
                fluorescence = np.random.normal(100, 10)
            else:  # Row C = test samples
                well_type = "sample"
                fluorescence = np.random.normal(500, 75)

            well = fp.Well(
                position=well_pos,
                fluorescence=fluorescence,
                well_type=well_type
            )
            plate.add_well(well)

    print(f"✅ Created plate with {len(plate.wells)} wells")

    # ======================================================================
    # 2. Use ANALYSIS FUNCTIONS (separated from data structure)
    # ======================================================================

    # Define control wells
    pos_controls = ['A1', 'A2', 'A3', 'A4']
    neg_controls = ['B1', 'B2', 'B3', 'B4']
    test_wells = ['C1', 'C2', 'C3', 'C4']

    print("\n📊 Statistical Analysis:")
    print("-" * 25)

    # Calculate CV for controls
    pos_cv = fp.analysis.calculate_cv(plate, pos_controls)
    neg_cv = fp.analysis.calculate_cv(plate, neg_controls)
    print(f"Positive Control CV: {pos_cv:.1f}%")
    print(f"Negative Control CV: {neg_cv:.1f}%")

    # Calculate Z-factor
    z_factor = fp.analysis.calculate_z_factor(plate, pos_controls, neg_controls)
    print(f"Z-factor: {z_factor:.3f}")

    # Detect outliers
    outliers = fp.analysis.detect_outliers(plate, test_wells, method='iqr')
    print(f"Outlier wells: {outliers}")

    print("\n🎯 Normalization Analysis:")
    print("-" * 27)

    # Normalize test wells to controls
    normalized = fp.analysis.normalize_to_controls(
        plate, test_wells, pos_controls, neg_controls
    )
    for well, norm_val in normalized.items():
        print(f"Well {well}: {norm_val:.1f}% of control")

    # Calculate percent inhibition
    inhibition = fp.analysis.percent_inhibition(plate, test_wells, pos_controls)
    for well, inhib_val in inhibition.items():
        print(f"Well {well}: {inhib_val:.1f}% inhibition")

    print("\n🔍 Quality Control:")
    print("-" * 19)

    # Generate QC report
    qc_report = fp.analysis.generate_qc_report(plate, pos_controls, neg_controls)
    print(f"QC Passed: {qc_report['qc_passed']}")
    print(f"Plate Mean Signal: {qc_report['plate_statistics']['mean_signal']:.0f}")
    print(f"Signal Ratio (Pos/Neg): {qc_report['control_validation']['signal_ratio']:.1f}")

    # ======================================================================
    # 3. Show the benefits of this architecture
    # ======================================================================

    print("\n🏗️ Architecture Benefits:")
    print("-" * 25)
    print("✅ Plate class stays focused on data storage")
    print("✅ Analysis functions are reusable and testable")
    print("✅ Easy to add new analysis without changing Plate")
    print("✅ Users import only the analysis they need")
    print("✅ Clear separation of concerns")

    # Example: Adding new analysis is easy
    print("\n📈 Easy to extend:")
    print("from fluoropy.analysis.statistics import calculate_cv")
    print("from fluoropy.analysis.normalization import z_score_normalize")
    print("from fluoropy.analysis.quality_control import check_edge_effects")

if __name__ == "__main__":
    main()
