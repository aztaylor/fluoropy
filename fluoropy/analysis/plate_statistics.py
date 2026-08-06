"""
Plate-level statistics: per-group summaries, outliers and z-scores.

These operate across a whole plate at one timepoint, answering "how do the
wells compare to each other" rather than "what does this sample do over time"
(which is Sample and SampleFrame territory).

They were methods on Plate. They only ever iterated wells and read
``well.time_series``, so they live here instead, following the same split
already used for plotting. Plate keeps thin delegating methods, so
``plate.calculate_zscore_normalization(...)`` continues to work unchanged --
these are equivalent:

    plate.calculate_zscore_normalization("OD600", timepoint_idx=10)
    fluoropy.analysis.calculate_zscore_normalization(plate, "OD600", 10)

For a z-score that is not dragged around by the outliers you are hunting, see
:func:`fluoropy.analysis.robust_z_score`.
"""

import logging
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from ..core.well import parse_well_id

logger = logging.getLogger(__name__)


def calculate_timepoint_statistics(plate, measurement_type: str, timepoint_idx: int,
                                 sample_types: Optional[List[str]] = None,
                                 exclude_blanks: bool = True,
                                 exclude_controls: bool = False) -> Dict[str, Dict[str, Any]]:
    """
    Calculate summary statistics for each sample type and concentration at a given timepoint.

    Parameters
    ----------
    measurement_type : str
        Type of measurement to analyze
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    sample_types : List[str], optional
        Specific sample types to analyze. If None, analyzes all sample types.
    exclude_blanks : bool, default True
        Whether to exclude blank wells from analysis
    exclude_controls : bool, default False
        Whether to exclude control wells from analysis

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary with 'sample_type_concentration' keys, and statistics dictionaries as values.
        Each statistics dictionary contains: 'mean', 'std', 'sem', 'count', 'min', 'max', 'median', 'q25', 'q75', 'iqr', 'outlier_wells'

    Examples
    --------
    >>> stats = plate.calculate_timepoint_statistics('OD600', timepoint_idx=10)
    >>> print(stats['sample_1_10.0']['mean'])  # Mean OD600 for sample_1 at 10.0 concentration
    >>> print(stats['sample_1_5.0']['outlier_wells'])  # List of outlier well IDs for sample_1 at 5.0 concentration
    """
    from collections import defaultdict

    # Group wells by (sample_type, concentration) with well information
    sample_conc_groups = defaultdict(list)

    for well in plate.wells.values():
        # Skip excluded wells
        if well.is_excluded():
            continue

        # Skip wells without the measurement
        if not (hasattr(well, 'time_series') and measurement_type in well.time_series):
            continue

        # Skip wells without enough timepoints
        time_series = well.time_series[measurement_type]
        if len(time_series) <= timepoint_idx:
            continue

        # Apply exclusion criteria
        if exclude_blanks and hasattr(well, 'is_blank') and well.is_blank:
            continue
        if exclude_controls and hasattr(well, 'is_control') and well.is_control:
            continue

        # Get sample type and concentration
        sample_type = getattr(well, 'sample_type', 'Unknown')
        concentration = getattr(well, 'concentration', 0.0)

        # Filter by specific sample types if provided
        if sample_types is not None and sample_type not in sample_types:
            continue

        # Create group key (sample_type, concentration)
        group_key = f"{sample_type}_{concentration}"

        # Add value and well information to group
        value = time_series[timepoint_idx]
        well_id = getattr(well, 'well_id', getattr(well, 'position', 'Unknown'))
        sample_conc_groups[group_key].append({'value': value, 'well_id': well_id, 'well': well})

    # Calculate statistics for each (sample_type, concentration) group
    statistics = {}

    for group_key, well_data_list in sample_conc_groups.items():
        if not well_data_list:
            continue

        # Extract values and well information
        values = [item['value'] for item in well_data_list]
        well_ids = [item['well_id'] for item in well_data_list]
        wells = [item['well'] for item in well_data_list]

        values_array = np.array(values)

        # Calculate basic statistics
        stats = {
            'mean': np.mean(values_array),
            'std': np.std(values_array, ddof=1) if len(values_array) > 1 else 0.0,
            'sem': np.std(values_array, ddof=1) / np.sqrt(len(values_array)) if len(values_array) > 1 else 0.0,
            'count': len(values_array),
            'min': np.min(values_array),
            'max': np.max(values_array),
            'median': np.median(values_array)
        }

        # Add quartiles and IQR
        stats['q25'] = np.percentile(values_array, 25)
        stats['q75'] = np.percentile(values_array, 75)
        stats['iqr'] = stats['q75'] - stats['q25']

        # Identify outliers using IQR method (values outside Q1 - 1.5*IQR or Q3 + 1.5*IQR)
        outlier_wells = []
        if len(values_array) > 2 and stats['iqr'] > 0:  # Need at least 3 values and non-zero IQR
            lower_bound = stats['q25'] - 1.5 * stats['iqr']
            upper_bound = stats['q75'] + 1.5 * stats['iqr']

            for i, value in enumerate(values):
                if value < lower_bound or value > upper_bound:
                    outlier_wells.append({
                        'well_id': well_ids[i],
                        'value': value,
                        'z_score': (value - stats['mean']) / stats['std'] if stats['std'] > 0 else 0
                    })

        stats['outlier_wells'] = outlier_wells
        stats['outlier_count'] = len(outlier_wells)

        statistics[group_key] = stats

    return statistics

def get_timepoint_summary_table(plate, measurement_type: str, timepoint_idx: int,
                              sample_types: Optional[List[str]] = None,
                              exclude_blanks: bool = True,
                              exclude_controls: bool = False,
                              include_outliers: bool = True) -> 'pd.DataFrame':
    """
    Get summary statistics as a formatted pandas DataFrame.

    Parameters
    ----------
    measurement_type : str
        Type of measurement to analyze
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    sample_types : List[str], optional
        Specific sample types to analyze
    exclude_blanks : bool, default True
        Whether to exclude blank wells
    exclude_controls : bool, default False
        Whether to exclude control wells
    include_outliers : bool, default True
        Whether to include outlier information in the DataFrame

    Returns
    -------
    pd.DataFrame
        DataFrame with sample types as index and statistics as columns.
        If include_outliers=True, includes 'outlier_count' and 'outlier_wells' columns.
    """
    stats = calculate_timepoint_statistics(
        plate, measurement_type, timepoint_idx, sample_types,
        exclude_blanks, exclude_controls
    )

    if not stats:
        return pd.DataFrame()

    # Convert to DataFrame
    df = pd.DataFrame.from_dict(stats, orient='index')

    # Round numeric columns
    numeric_columns = ['mean', 'std', 'sem', 'min', 'max', 'median', 'q25', 'q75', 'iqr']
    for col in numeric_columns:
        if col in df.columns:
            df[col] = df[col].round(4)

    if include_outliers:
        # Format outlier wells for better display
        if 'outlier_wells' in df.columns:
            df['outlier_wells_formatted'] = df['outlier_wells'].apply(
                lambda x: '; '.join([f"{well['well_id']}({well['value']:.2f})" for well in x]) if x else 'None'
            )

            # Keep outlier_count for quick reference
            df['outlier_count'] = df['outlier_count'] if 'outlier_count' in df.columns else 0

            # Reorder columns to put outlier info at the end
            base_cols = [col for col in df.columns if col not in ['outlier_wells', 'outlier_wells_formatted', 'outlier_count']]
            outlier_cols = ['outlier_count', 'outlier_wells_formatted']
            df = df[base_cols + outlier_cols]
    else:
        # Remove outlier columns if not wanted
        outlier_cols = ['outlier_wells', 'outlier_count', 'outlier_wells_formatted']
        df = df.drop(columns=[col for col in outlier_cols if col in df.columns])

    # Sort by sample type
    df = df.sort_index()

    return df

def get_outlier_wells(plate, measurement_type: str, timepoint_idx: int,
                     sample_types: Optional[List[str]] = None,
                     exclude_blanks: bool = True,
                     exclude_controls: bool = False) -> Dict[str, List[Dict[str, Any]]]:
    """
    Get outlier wells for each sample type and concentration at a given timepoint.

    Parameters
    ----------
    measurement_type : str
        Type of measurement to analyze
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    sample_types : List[str], optional
        Specific sample types to analyze
    exclude_blanks : bool, default True
        Whether to exclude blank wells
    exclude_controls : bool, default False
        Whether to exclude control wells

    Returns
    -------
    Dict[str, List[Dict[str, Any]]]
        Dictionary with 'sample_type_concentration' keys and lists of outlier well information.
        Each outlier dictionary contains: 'well_id', 'value', 'z_score'

    Examples
    --------
    >>> outliers = plate.get_outlier_wells('OD600', timepoint_idx=10)
    >>> print(outliers['sample_1_10.0'])  # List of outlier wells for sample_1 at 10.0 concentration
    """
    stats = calculate_timepoint_statistics(
        plate, measurement_type, timepoint_idx, sample_types,
        exclude_blanks, exclude_controls
    )

    outliers = {}
    for group_key, group_stats in stats.items():
        if 'outlier_wells' in group_stats:
            outliers[group_key] = group_stats['outlier_wells']

    return outliers

def calculate_zscore_normalization(plate, measurement_type: str, timepoint_idx: int,
                                 exclude_blanks: bool = True,
                                 exclude_controls: bool = False) -> Dict[str, float]:
    """
    Calculate z-score normalization for all wells at a specific timepoint.

    Z-score = (value - mean) / std_dev

    This normalizes values across the entire plate, making it easy to identify
    wells that deviate significantly from the plate average.

    Parameters
    ----------
    measurement_type : str
        Type of measurement to normalize
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    exclude_blanks : bool, default True
        Whether to exclude blank wells from the calculation of plate statistics
    exclude_controls : bool, default False
        Whether to exclude control wells from the calculation of plate statistics

    Returns
    -------
    Dict[str, float]
        Dictionary with well IDs as keys and z-scores as values.
        Wells excluded from calculation will not be included in the result.

    Examples
    --------
    >>> z_scores = plate.calculate_zscore_normalization('OD600', timepoint_idx=10)
    >>> print(z_scores['A1'])  # Z-score for well A1
    >>> extreme_wells = {k: v for k, v in z_scores.items() if abs(v) > 2}  # Wells with |z| > 2
    """
    # Collect all values for plate-wide statistics
    all_values = []
    well_values = {}

    for well in plate.wells.values():
        # Skip excluded wells
        if well.is_excluded():
            continue

        # Skip wells without the measurement
        if not (hasattr(well, 'time_series') and measurement_type in well.time_series):
            continue

        # Skip wells without enough timepoints
        time_series = well.time_series[measurement_type]
        if len(time_series) <= timepoint_idx:
            continue

        # Apply exclusion criteria for plate statistics calculation
        if exclude_blanks and hasattr(well, 'is_blank') and well.is_blank:
            continue
        if exclude_controls and hasattr(well, 'is_control') and well.is_control:
            continue

        # Get well ID and value
        well_id = getattr(well, 'well_id', getattr(well, 'position', 'Unknown'))
        value = time_series[timepoint_idx]

        all_values.append(value)
        well_values[well_id] = value

    if len(all_values) < 2:
        # Need at least 2 values to calculate standard deviation
        return {}

    # Calculate plate-wide statistics
    plate_mean = np.mean(all_values)
    plate_std = np.std(all_values, ddof=1)

    if plate_std == 0:
        # If standard deviation is 0, all values are the same
        return {well_id: 0.0 for well_id in well_values.keys()}

    # Calculate z-scores
    z_scores = {}
    for well_id, value in well_values.items():
        z_scores[well_id] = (value - plate_mean) / plate_std

    return z_scores

def apply_zscore_normalization(plate, measurement_type: str, timepoint_idx: int,
                             exclude_blanks: bool = True,
                             exclude_controls: bool = False,
                             store_in_metadata: bool = True) -> Dict[str, float]:
    """
    Apply z-score normalization and optionally store results in well metadata.

    Parameters
    ----------
    measurement_type : str
        Type of measurement to normalize
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    exclude_blanks : bool, default True
        Whether to exclude blank wells from plate statistics calculation
    exclude_controls : bool, default False
        Whether to exclude control wells from plate statistics calculation
    store_in_metadata : bool, default True
        Whether to store z-scores in well metadata for later access

    Returns
    -------
    Dict[str, float]
        Dictionary with well IDs as keys and z-scores as values

    Examples
    --------
    >>> z_scores = plate.apply_zscore_normalization('OD600', timepoint_idx=10)
    >>> # Z-scores are now stored in each well's metadata
    >>> well_a1 = plate['A1']
    >>> print(well_a1.metadata.get('zscore_OD600_tp10'))  # Access stored z-score
    """
    z_scores = calculate_zscore_normalization(
        plate, measurement_type, timepoint_idx, exclude_blanks, exclude_controls
    )

    if store_in_metadata:
        # Store z-scores in well metadata
        metadata_key = f"zscore_{measurement_type}_tp{timepoint_idx}"

        for well_id, z_score in z_scores.items():
            well = plate.wells.get(well_id)
            if well:
                if not hasattr(well, 'metadata'):
                    well.metadata = {}
                well.metadata[metadata_key] = z_score

    return z_scores

def get_zscore_matrix(plate, measurement_type: str, timepoint_idx: int,
                     exclude_blanks: bool = True,
                     exclude_controls: bool = False) -> np.ndarray:
    """
    Get z-scores as a 2D matrix matching the plate layout for visualization.

    Parameters
    ----------
    measurement_type : str
        Type of measurement to normalize
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    exclude_blanks : bool, default True
        Whether to exclude blank wells from plate statistics calculation
    exclude_controls : bool, default False
        Whether to exclude control wells from plate statistics calculation

    Returns
    -------
    np.ndarray
        2D array of z-scores with shape (rows, cols) matching plate layout.
        Wells excluded from analysis will have NaN values.

    Examples
    --------
    >>> z_matrix = plate.get_zscore_matrix('OD600', timepoint_idx=10)
    >>> import matplotlib.pyplot as plt
    >>> plt.imshow(z_matrix, cmap='RdBu_r', vmin=-3, vmax=3)
    >>> plt.colorbar(label='Z-score')
    >>> plt.title('Plate Z-score Heatmap')
    """
    z_scores = calculate_zscore_normalization(
        plate, measurement_type, timepoint_idx, exclude_blanks, exclude_controls
    )

    # Create matrix with NaN for missing values
    z_matrix = np.full((plate.rows, plate.cols), np.nan)

    for well_id, z_score in z_scores.items():
        well = plate.wells.get(well_id)
        if well:
            # Parse row and column from well_id
            row_idx, col_idx = parse_well_id(well_id)

            if 0 <= row_idx < plate.rows and 0 <= col_idx < plate.cols:
                z_matrix[row_idx, col_idx] = z_score

    return z_matrix

