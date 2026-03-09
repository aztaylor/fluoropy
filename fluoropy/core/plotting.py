"""
Plotting functions for fluoropy core objects.

This module contains all visualization code extracted from Plate and SampleFrame
classes. Functions accept core objects as their first argument.

Requires matplotlib (optional dependency, install with: pip install fluoropy[viz])
"""

from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import matplotlib.figure


def _require_matplotlib():
    """Import and return matplotlib.pyplot, raising a clear error if unavailable."""
    try:
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install fluoropy[viz]"
        )


# ==========================================================================
# Plate plotting functions
# ==========================================================================

def plot_timeseries_grid(plate, measurement_type: str, figsize: tuple = (15, 10),
                         title: Optional[str] = None,
                         show_sample_info: bool = True) -> 'matplotlib.figure.Figure':
    """
    Plot raw timeseries data for each well in an 8x12 grid matching plate layout.

    Parameters
    ----------
    plate : Plate
        Plate object to plot
    measurement_type : str
        Type of measurement to plot (e.g., 'OD600', 'GFP')
    figsize : tuple, default (15, 10)
        Figure size in inches (width, height)
    title : str, optional
        Overall figure title. If None, uses measurement type and plate name.
    show_sample_info : bool, default True
        Whether to show sample type and concentration in subplot titles

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object
    """
    plt = _require_matplotlib()

    fig, axes = plt.subplots(8, 12, figsize=figsize, sharex=True, sharey=True)
    fig.suptitle(title or f"{measurement_type} Timeseries - {plate.name or 'Plate'}",
                fontsize=14, fontweight='bold')

    for row in range(8):
        for col in range(12):
            ax = axes[row, col]
            well_id = f"{chr(ord('A') + row)}{col + 1}"
            well = plate.wells.get(well_id)

            if well and hasattr(well, 'time_series') and measurement_type in well.time_series:
                y_data = well.time_series[measurement_type]
                x_data = well.time_points if hasattr(well, 'time_points') and well.time_points is not None else range(len(y_data))

                ax.plot(x_data, y_data, 'b-', linewidth=1.5, alpha=0.8)

                if hasattr(well, 'is_blank') and well.is_blank:
                    ax.set_facecolor('#f0f0f0')
                elif hasattr(well, 'is_control') and well.is_control:
                    ax.set_facecolor('#fff5f5')
                else:
                    ax.set_facecolor('white')

                if show_sample_info and hasattr(well, 'sample_type'):
                    sample_info = well.sample_type or "Unknown"
                    if hasattr(well, 'concentration') and well.concentration is not None:
                        sample_info += f"\n{well.concentration}"
                    ax.set_title(f"{well_id}\n{sample_info}", fontsize=8, pad=2)
                else:
                    ax.set_title(well_id, fontsize=8, pad=2)
            else:
                ax.set_facecolor('#e0e0e0')
                ax.set_title(well_id, fontsize=8, pad=2)
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=8, alpha=0.6)

            ax.tick_params(labelsize=6)
            ax.grid(True, alpha=0.3)

    fig.text(0.5, 0.02, 'Time', ha='center', fontsize=12)
    fig.text(0.02, 0.5, measurement_type, va='center', rotation='vertical', fontsize=12)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.07, left=0.05, right=0.98)

    return fig


def plot_zscore_heatmap(plate, measurement_type: str, timepoint_idx: int,
                        exclude_blanks: bool = True,
                        exclude_controls: bool = False,
                        figsize: tuple = (10, 6),
                        title: Optional[str] = None,
                        vmin: float = -3,
                        vmax: float = 3) -> 'matplotlib.figure.Figure':
    """
    Plot a heatmap of z-scores across the plate.

    Parameters
    ----------
    plate : Plate
        Plate object to plot
    measurement_type : str
        Type of measurement to plot
    timepoint_idx : int
        Index of the timepoint to analyze (0-based)
    exclude_blanks : bool, default True
        Whether to exclude blank wells from plate statistics calculation
    exclude_controls : bool, default False
        Whether to exclude control wells from plate statistics calculation
    figsize : tuple, default (10, 6)
        Figure size in inches (width, height)
    title : str, optional
        Plot title. If None, generates automatic title.
    vmin : float, default -3
        Minimum value for color scale
    vmax : float, default 3
        Maximum value for color scale

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object
    """
    plt = _require_matplotlib()

    z_matrix = plate.get_zscore_matrix(
        measurement_type, timepoint_idx, exclude_blanks, exclude_controls
    )

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(z_matrix, cmap='RdBu_r', vmin=vmin, vmax=vmax, aspect='auto')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Z-score', rotation=270, labelpad=20)

    ax.set_xlabel('Column')
    ax.set_ylabel('Row')

    ax.set_xticks(range(plate.cols))
    ax.set_xticklabels([str(i+1) for i in range(plate.cols)])
    ax.set_yticks(range(plate.rows))
    ax.set_yticklabels([chr(ord('A') + i) for i in range(plate.rows)])

    if title is None:
        title = f'Z-score Heatmap: {measurement_type} (Timepoint {timepoint_idx})'
    ax.set_title(title)

    for i in range(plate.rows):
        for j in range(plate.cols):
            z_val = z_matrix[i, j]
            if not np.isnan(z_val) and abs(z_val) > 2:
                ax.text(j, i, f'{z_val:.1f}', ha='center', va='center',
                       color='white' if abs(z_val) > 2.5 else 'black',
                       fontweight='bold', fontsize=8)

    plt.tight_layout()
    return fig


# ==========================================================================
# SampleFrame plotting functions
# ==========================================================================

def plot_fold_change_dose_response(frame, timepoint_idx: int,
                                   sample_ids: Optional[List[str]] = None,
                                   x_logscale: bool = True,
                                   y_logscale: bool = True,
                                   figsize: Optional[Tuple[int, int]] = None,
                                   title: Optional[str] = None,
                                   ylabel: str = "Fold Change",
                                   xlabel: str = "Concentration",
                                   concentration_idx_range: Optional[Tuple[int, int]] = None):
    """
    Plot dose-response curve: fold change vs concentration at a specific timepoint.

    Parameters
    ----------
    frame : SampleFrame
        SampleFrame object containing fold change data
    timepoint_idx : int
        Index of the timepoint to plot (0-based)
    sample_ids : List[str], optional
        List of sample IDs to plot. If None, plots all test samples.
    x_logscale : bool, default True
        If True, use log scale for x-axis (concentration)
    y_logscale : bool, default True
        If True, use log scale for y-axis (fold change)
    figsize : Tuple[int, int], optional
        Figure size (width, height) in inches
    title : str, optional
        Figure title. If None, auto-generates based on timepoint index
    ylabel : str, default "Fold Change"
        Y-axis label
    xlabel : str, default "Concentration"
        X-axis label
    concentration_idx_range : Tuple[int, int], optional
        Tuple of (start_idx, end_idx) to plot only a range of concentrations.

    Returns
    -------
    Tuple[Figure, Axes]
        Figure and axes objects
    """
    plt = _require_matplotlib()

    if not sample_ids:
        sample_ids = frame.get_test_samples()

    if not sample_ids:
        raise ValueError("No test samples found")

    for sid in sample_ids:
        if sid not in frame.samples:
            raise ValueError(f"Sample ID '{sid}' not found")

    if figsize is None:
        figsize = (10, 6)
    fig, ax = plt.subplots(figsize=figsize)

    if title is None:
        title = f"Dose-Response Curve at Timepoint {timepoint_idx}"
    ax.set_title(title, fontsize=14, fontweight='bold')

    colors = plt.cm.Set1(np.linspace(0, 1, len(sample_ids)))

    for color_idx, sample_id in enumerate(sample_ids):
        sample = frame.samples[sample_id]

        if not hasattr(sample, 'fold_change_mean'):
            print(f"Warning: No fold_change_mean data for {sample_id}")
            continue

        mean_dict = sample.fold_change_mean
        error_dict = sample.fold_change_error

        all_concentrations = sorted([c for c in mean_dict.keys() if c != 0.0])

        if concentration_idx_range is not None:
            start_idx, end_idx = concentration_idx_range
            concentrations = all_concentrations[start_idx:end_idx+1]
        else:
            concentrations = all_concentrations

        if not concentrations:
            print(f"Warning: No concentrations in specified range for {sample_id}")
            continue

        fold_changes = []
        errors = []

        for conc in concentrations:
            mean_array = mean_dict[conc]
            error_array = error_dict[conc]

            if timepoint_idx >= len(mean_array):
                raise IndexError(f"Timepoint index {timepoint_idx} out of range (max {len(mean_array)-1})")

            fold_changes.append(mean_array[timepoint_idx])
            errors.append(error_array[timepoint_idx])

        ax.errorbar(int(concentrations), fold_changes, yerr=errors,
                   marker='o', linestyle='-', linewidth=2, markersize=8,
                   capsize=5, capthick=2, label=sample_id, color=colors[color_idx])

    if x_logscale:
        ax.set_xscale('log')
    if y_logscale:
        ax.set_yscale('log')

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best', fontsize=11)

    plt.tight_layout()
    return fig, ax


def plot_replicate_time_series(frame,
                               measurement: str,
                               sample_ids: Optional[List[str]] = None,
                               show_mean: bool = True,
                               figsize: Optional[Tuple[int, int]] = None,
                               title: Optional[str] = None,
                               ylabel: Optional[str] = None,
                               xlabel: str = "Time (hours)",
                               sharey: bool = True,
                               max_cols: int = 4) -> Tuple['matplotlib.figure.Figure', Dict[str, 'matplotlib.axes.Axes']]:
    """
    Plot time series curves for replicates of samples at each concentration.

    Creates a subplot for each sample-concentration combination, displaying the
    individual replicate curves (wells) within that subplot.

    Parameters
    ----------
    frame : SampleFrame
        SampleFrame object containing samples with well data
    measurement : str
        Measurement type to plot (e.g., 'OD600', 'GFP')
    sample_ids : List[str], optional
        List of sample IDs to include. If None, includes all test samples.
    show_mean : bool, default True
        Whether to overlay the mean curve across replicates with error band
    figsize : Tuple[int, int], optional
        Figure size (width, height) in inches
    title : str, optional
        Figure title
    ylabel : str, optional
        Y-axis label
    xlabel : str, default "Time (hours)"
        X-axis label
    sharey : bool, default True
        Whether subplots share y-axis

    Returns
    -------
    Tuple[Figure, Dict[str, Axes]]
        Figure object and dictionary mapping subplot keys to Axes objects
    """
    plt = _require_matplotlib()

    if not sample_ids:
        sample_ids = frame.get_test_samples()

    if not sample_ids:
        raise ValueError("No test samples found. Check your sample data.")

    for sid in sample_ids:
        if sid not in frame.samples:
            raise ValueError(f"Sample ID '{sid}' not found in SampleFrame")

    # Check first that we have valid wells with the measurement
    measurement_found = False
    wells_with_data = 0

    for sample_id in sample_ids:
        sample = frame.samples[sample_id]
        if not sample.wells:
            raise RuntimeError(f"Sample '{sample_id}' has no wells")

        for well in sample.wells:
            if measurement in well.time_series and not well.is_excluded():
                measurement_found = True
                if well.time_points is not None:
                    wells_with_data += 1

    if not measurement_found:
        available_measurements = set()
        for sample_id in sample_ids:
            for well in frame.samples[sample_id].wells:
                available_measurements.update(well.time_series.keys())
        raise ValueError(f"Measurement '{measurement}' not found in wells. Available: {list(available_measurements) if available_measurements else 'none'}")

    if wells_with_data == 0:
        raise RuntimeError("No wells with valid time_points data found")

    # Collect all (sample_id, concentration) pairs with their wells
    subplot_data = {}

    for sample_id in sample_ids:
        sample = frame.samples[sample_id]
        conc_groups = {}
        for well in sample.wells:
            if not well.is_excluded() and measurement in well.time_series:
                conc = well.concentration
                if conc not in conc_groups:
                    conc_groups[conc] = []
                conc_groups[conc].append(well)

        for conc, wells in conc_groups.items():
            if wells:
                seen_wells = {}
                unique_wells = []
                for well in wells:
                    well_key = id(well)
                    if well_key not in seen_wells:
                        seen_wells[well_key] = True
                        unique_wells.append(well)

                if unique_wells:
                    key = (sample_id, conc)
                    subplot_data[key] = unique_wells

    if not subplot_data:
        raise ValueError("No wells found with valid data")

    sorted_keys = sorted(subplot_data.keys(), key=lambda x: (x[0], x[1]))
    n_subplots = len(sorted_keys)

    if figsize is None:
        n_cols = int(np.ceil(np.sqrt(n_subplots / 1.4)))
        n_cols = max(max_cols, min(n_cols, 8))
        n_rows = int(np.ceil(n_subplots / n_cols))
        figsize = (n_cols * 3.5, n_rows * 3.5)
    else:
        n_cols = int(np.ceil(np.sqrt(n_subplots / 1.4)))
        n_cols = max(max_cols, min(n_cols, 8))
        n_rows = int(np.ceil(n_subplots / n_cols))

    fig, axes_array = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False, sharey=sharey)
    axes_flat = axes_array.flatten()

    if title is None:
        title = f"Replicate Time Series - {measurement}"
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)

    axes_dict = {}

    for subplot_idx, (sample_id, conc) in enumerate(sorted_keys):
        ax = axes_flat[subplot_idx]
        wells = subplot_data[(sample_id, conc)]

        for well in wells:
            if measurement in well.time_series:
                replicate_data = well.time_series[measurement]
                if well.time_points is not None:
                    time = well.time_points
                    # Trim trailing zeros from pre-allocated import arrays
                    nonzero = np.nonzero(time)[0]
                    n = int(nonzero[-1]) + 1 if len(nonzero) else len(time)
                    time = time[:n]
                    replicate_data = replicate_data[:n]
                else:
                    time = np.arange(len(replicate_data))
                label = getattr(well, 'plate_id', well.well_id)
                ax.plot(time, replicate_data, '-', alpha=0.6,
                       linewidth=1.5, label=label, zorder=2)

        if show_mean:
            replicate_arrays = np.column_stack([w.time_series[measurement] for w in wells])
            mean_data = np.mean(replicate_arrays, axis=1)
            error_data = np.std(replicate_arrays, axis=1)
            ax.plot(time, mean_data, 'o-', linewidth=2.5,
                   markersize=6, label=f"{sample_id} mean", zorder=3)

            ax.fill_between(time,
                           mean_data - error_data,
                           mean_data + error_data,
                           alpha=0.15, zorder=1)

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel if ylabel else measurement, fontsize=10)
        ax.set_title(f"{sample_id} [{conc}]", fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='best', fontsize=8, title="Replicates", title_fontsize=9)

        subplot_key = f"({sample_id}, {conc})"
        axes_dict[subplot_key] = ax

    for idx in range(n_subplots, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.tight_layout()
    return fig, axes_dict


def plot_dose_response_with_hill_fit(frame, timepoint_idx: int,
                                     sample_ids: Optional[List[str]] = None,
                                     x_logscale: bool = True,
                                     y_logscale: bool = True,
                                     figsize: Optional[Tuple[int, int]] = None,
                                     title: Optional[str] = None,
                                     ylabel: str = "Fold Change",
                                     xlabel: str = "Concentration",
                                     concentration_idx_range: Optional[Tuple[int, int]] = None,
                                     num_points: int = 100):
    """
    Plot dose-response curve with Hill function fit overlay.

    Parameters
    ----------
    frame : SampleFrame
        SampleFrame object containing fold change data
    timepoint_idx : int
        Index of the timepoint to plot (0-based)
    sample_ids : List[str], optional
        List of sample IDs to plot. If None, plots all test samples.
    x_logscale : bool, default True
        If True, use log scale for x-axis
    y_logscale : bool, default True
        If True, use log scale for y-axis
    figsize : Tuple[int, int], optional
        Figure size (width, height) in inches
    title : str, optional
        Figure title
    ylabel : str, default "Fold Change"
        Y-axis label
    xlabel : str, default "Concentration"
        X-axis label
    concentration_idx_range : Tuple[int, int], optional
        Tuple of (start_idx, end_idx) to plot only a range of concentrations.
    num_points : int, default 100
        Number of points for the fitted curve

    Returns
    -------
    Tuple[Figure, Axes]
        Figure and axes objects
    """
    plt = _require_matplotlib()

    fit_results = frame.calculate_hill_fits(timepoint_idx, sample_ids, concentration_idx_range)

    if not sample_ids:
        sample_ids = frame.get_test_samples()

    if not sample_ids:
        raise ValueError("No test samples found")

    if figsize is None:
        figsize = (12, 7)
    fig, ax = plt.subplots(figsize=figsize)

    if title is None:
        title = f"Dose-Response Curve (Hill Fit) at Timepoint {timepoint_idx}"
    ax.set_title(title, fontsize=14, fontweight='bold')

    def hill_function(conc, ec50, hill, min_val, max_val):
        return min_val + (max_val - min_val) * (conc ** hill) / (ec50 ** hill + conc ** hill)

    colors = plt.cm.Set1(np.linspace(0, 1, len(sample_ids)))

    for color_idx, sample_id in enumerate(sample_ids):
        sample = frame.samples[sample_id]

        if not hasattr(sample, 'fold_change_mean'):
            print(f"Warning: No fold_change_mean data for {sample_id}")
            continue

        mean_dict = sample.fold_change_mean
        error_dict = sample.fold_change_error

        all_concentrations = sorted([c for c in mean_dict.keys() if c != 0.0])

        if concentration_idx_range is not None:
            start_idx, end_idx = concentration_idx_range
            concentrations = all_concentrations[start_idx:end_idx+1]
        else:
            concentrations = all_concentrations

        if not concentrations:
            continue

        fold_changes = []
        errors = []
        for conc in concentrations:
            mean_array = mean_dict[conc]
            error_array = error_dict[conc]
            if timepoint_idx >= len(mean_array):
                raise IndexError(f"Timepoint index {timepoint_idx} out of range")
            fold_changes.append(mean_array[timepoint_idx])
            errors.append(error_array[timepoint_idx])

        ax.errorbar(concentrations, fold_changes, yerr=errors,
                   marker='o', linestyle='', linewidth=2, markersize=8,
                   capsize=5, capthick=2, label=sample_id, color=colors[color_idx],
                   alpha=0.7)

        if sample_id in fit_results:
            params = fit_results[sample_id]
            ec50 = params['ec50']
            hill = params['hill']
            min_val = params['min']
            max_val = params['max']
            r_sq = params['r_squared']

            conc_min, conc_max = np.min(concentrations), np.max(concentrations)
            if x_logscale:
                conc_smooth = np.logspace(np.log10(conc_min), np.log10(conc_max), num_points)
            else:
                conc_smooth = np.linspace(conc_min, conc_max, num_points)

            fc_smooth = hill_function(conc_smooth, ec50, hill, min_val, max_val)

            ax.plot(conc_smooth, fc_smooth, '-', linewidth=2.5,
                   color=colors[color_idx], alpha=0.9,
                   label=f"{sample_id} fit (EC50={ec50:.2e}, n={hill:.2f}, R²={r_sq:.3f})")

    if x_logscale:
        ax.set_xscale('log')
    if y_logscale:
        ax.set_yscale('log')

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best', fontsize=10)

    plt.tight_layout()
    return fig, ax


def plot_mean_normalized_data(frame, measurement: str,
                              sample_ids: Optional[List[str]] = None,
                              figsize: Optional[Tuple[int, int]] = None,
                              title: Optional[str] = None,
                              legend_title: str = "Concentration",
                              ylabel: Optional[str] = None,
                              xlabel: str = "Time (hours)",
                              show_error: bool = True,
                              error_type: str = 'std',
                              error_alpha: float = 0.1,
                              concentration_idx_range: Optional[Tuple[int, int]] = None):
    """
    Plot OD normalized fluorescence time series for each sample as subplots.

    Parameters
    ----------
    frame : SampleFrame
        SampleFrame object containing normalized data
    measurement : str
        Measurement type to plot (e.g., 'GFP', 'RFP')
    sample_ids : List[str], optional
        List of sample IDs to include. If None, includes all non-blank samples.
    figsize : Tuple[int, int], optional
        Figure size (width, height) in inches
    title : str, optional
        Figure title
    legend_title : str, default "Concentration"
        Title for the legend
    ylabel : str, optional
        Y-axis label
    xlabel : str, default "Time (hours)"
        X-axis label
    show_error : bool, default True
        Whether to show error bands around the mean curves
    error_type : str, default 'std'
        Type of error to display: 'std' or 'sem'
    error_alpha : float, default 0.1
        Transparency level for the error bands
    concentration_idx_range : Tuple[int, int], optional
        Tuple of (start_idx, end_idx) to plot only a range of concentrations.

    Returns
    -------
    Tuple[Figure, Dict[str, Axes]]
        Figure object and dictionary mapping sample IDs to Axes objects
    """
    plt = _require_matplotlib()

    if not sample_ids:
        sample_ids = [sid for sid in frame.samples.keys()
                     if not frame.samples[sid].is_blank]

    if not sample_ids:
        raise ValueError("No samples found (excluding blanks)")

    for sid in sample_ids:
        if sid not in frame.samples:
            raise ValueError(f"Sample ID '{sid}' not found in SampleFrame")

    has_normalized_data = False
    for sample_id in sample_ids:
        sample = frame.samples[sample_id]
        if hasattr(sample, 'normalized_data_mean') and sample.normalized_data_mean:
            if measurement in sample.normalized_data_mean:
                has_normalized_data = True
                break

    if not has_normalized_data:
        raise RuntimeError(f"No normalized timeseries data found for measurement '{measurement}'. "
                         f"Call calculate_normalized_timeseries() first.")

    n_samples = len(sample_ids)
    if figsize is None:
        n_cols = int(np.ceil(np.sqrt(n_samples)))
        n_cols = max(1, min(n_cols, 5))
        n_rows = int(np.ceil(n_samples / n_cols))
        figsize = (n_cols * 4, n_rows * 3.5)
    else:
        n_cols = int(np.ceil(np.sqrt(n_samples)))
        n_cols = max(1, min(n_cols, 5))
        n_rows = int(np.ceil(n_samples / n_cols))

    fig, axes_array = plt.subplots(n_rows, n_cols, figsize=figsize,
                                   squeeze=False, sharey=True)
    axes_flat = axes_array.flatten()

    if title is None:
        title = f"OD Normalized {measurement} Time Series"
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.995)

    axes_dict = {}

    colors = plt.cm.viridis(np.linspace(0, 1, 10))

    for subplot_idx, sample_id in enumerate(sample_ids):
        ax = axes_flat[subplot_idx]
        sample = frame.samples[sample_id]

        if not hasattr(sample, 'normalized_data_mean') or not sample.normalized_data_mean:
            ax.text(0.5, 0.5, f"No data for {sample_id}",
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_visible(True)
            axes_dict[sample_id] = ax
            continue

        if measurement not in sample.normalized_data_mean:
            ax.text(0.5, 0.5, f"Measurement {measurement}\nnot found",
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_visible(True)
            axes_dict[sample_id] = ax
            continue

        if hasattr(sample, 'normalized_data_mean') and sample.normalized_data_mean:
            print(f"Normalized statistics found for {sample_id} with shape {sample.normalized_data_mean.get(measurement).shape if sample.normalized_data_mean.get(measurement) is not None else 'N/A'}")
            mean_data = sample.normalized_data_mean.get(measurement)
            error_data = sample.normalized_data_error.get(measurement) if show_error else None
        else:
            print(f"Warning: No statistics found for {sample_id}, using raw normalized data")
            mean_data = sample.normalized_timeseries[measurement]
            error_data = None

        if mean_data is None or mean_data.size == 0:
            ax.text(0.5, 0.5, f"No normalized data\nfor {sample_id}",
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_visible(True)
            axes_dict[sample_id] = ax
            continue

        time = sample.time
        if time is None:
            time = np.arange(mean_data.shape[0])

        concentrations = sample.get_concentrations()

        if concentration_idx_range is not None:
            start_idx, end_idx = concentration_idx_range
            conc_indices = np.arange(len(concentrations))[start_idx:end_idx+1]
            concentrations = concentrations[conc_indices]
        else:
            conc_indices = np.arange(len(concentrations))

        for idx, (conc_idx, conc) in enumerate(zip(conc_indices, concentrations)):
            color = colors[idx % len(colors)]

            if len(mean_data.shape) == 2:
                conc_mean = mean_data[:, conc_idx]
                conc_error = error_data[:, conc_idx] if error_data is not None else None
            else:
                conc_mean = mean_data
                conc_error = error_data if error_data is not None else None

            ax.scatter(time,
                       conc_mean,
                       s=1,
                       color=color,
                       label=f"[{conc}]",
                       zorder=3)

            if show_error and conc_error is not None:
                ax.fill_between(time,
                               conc_mean - conc_error,
                               conc_mean + conc_error,
                               alpha=error_alpha, color=color, zorder=1)

        ax.set_xlabel(xlabel, fontsize=10)
        if ylabel is None:
            ax.set_ylabel(f"OD Normalized {measurement}", fontsize=10)
        else:
            ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(sample_id, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        if len(concentrations) > 0:
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels,
                      loc='center',
                      bbox_to_anchor=(1.15, 0.5),
                      fontsize=9,
                      title=legend_title,
                      title_fontsize=10)

        axes_dict[sample_id] = ax

    for idx in range(len(sample_ids), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.tight_layout()
    return fig, axes_dict
