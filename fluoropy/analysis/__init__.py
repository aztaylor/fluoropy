"""
Analysis functions that complement the core data classes.

Scope: this package holds only what Plate, Sample and SampleFrame do not
already provide. Replicate statistics, IQR outliers, plate-wide z-scores,
fold change and Hill fits all live on the core classes -- use those:

    Sample.calculate_statistics()             replicate mean/std/sem
    Plate.calculate_timepoint_statistics()    per-group stats, IQR outliers
    Plate.calculate_zscore_normalization()    plate-wide z-score
    SampleFrame.calculate_fold_change()       fold change vs control
    SampleFrame.calculate_hill_fits()         dose-response fitting

Every function here takes a Plate, a measurement key and a timepoint index
(default -1, the endpoint), and reads from ``well.time_series``.
"""

from .normalization import (
    normalize_to_controls,
    percent_inhibition,
    robust_z_score,
)
from .quality import (
    calculate_signal_to_noise,
    calculate_z_factor,
    check_edge_effects,
)

__all__ = [
    # Assay quality
    "calculate_z_factor",
    "calculate_signal_to_noise",
    "check_edge_effects",
    # Normalization
    "normalize_to_controls",
    "percent_inhibition",
    "robust_z_score",
]
