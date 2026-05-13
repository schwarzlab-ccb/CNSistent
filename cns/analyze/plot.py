import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, patches as mpatches
from collections.abc import Sequence
from itertools import accumulate

from cns.utils.cytobands import cytoband_color
from cns.utils.gaps import gap_color
from cns.utils import get_cn_cols, hg19


def _get_CN_color(cn, min_cn, max_cn):
    val_range = max_cn - min_cn
    if np.isnan(cn):
        return (0.5, 0.5, 0.5)
    if cn == 0:
        return (1, 0, 0)  # Red for zero
    else:
        ratio = 1 - min(1, (cn - min_cn) / val_range)
        return (ratio, ratio, 1)  # Blue gradient based on cn


def _get_CN_color_vector(min_cn, max_cn):
    # Vectorize with otypes for float output
    vectorized_func = np.vectorize(lambda cn: _get_CN_color(cn, min_cn, max_cn), otypes=[float, float, float])
    return lambda cn_array: np.stack(vectorized_func(cn_array), axis=-1)


def _get_start(chrom, assembly):
    return assembly.chr_starts[chrom]


def _get_start_vector(assembly):
    return np.vectorize(lambda chrom: _get_start(chrom, assembly), otypes=[np.uint32])


def plot_lines(ax, cns_df, cn_column, color="green", label=None, alpha=1.0, size=1, assembly=hg19, linestyle="-"):
    """
    Plots consecutive segments as lines on the given axis - centers of each segment are used as endpoints of each line.

    NOTE: A single segment will not be plotted, at least two segments must exist.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the lines.
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for copy number data.
    color : str, optional
        Color of the lines. Default is "green".
    label : str, optional
        Label for the lines. Default is None.
    alpha : float, optional
        Alpha value for the lines. Default is 1.
    size : float, optional
        Line width. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted lines.
    """
    f_start_pos = _get_start_vector(assembly)
    is_consecutive = cns_df["start"] - cns_df["end"].shift(1) != 0
    # plot consecutive segments
    for _, group_df in cns_df.groupby(is_consecutive.cumsum()):
        length = group_df["end"] - group_df["start"]
        x = group_df["start"] + length / 2 + f_start_pos(group_df["chrom"])
        ax.plot(x, group_df[cn_column], c=color, linewidth=size, label=label, alpha=alpha, linestyle=linestyle)
        label = None  # only use label for the first segment
    return ax


def plot_steps(ax, cns_df, cn_column, color="green", label=None, alpha=1.0, size=1, assembly=hg19, linestyle="-"):
    f_start_pos = _get_start_vector(assembly)
    is_consecutive = cns_df["start"] - cns_df["end"].shift(1) != 0
    # plot consecutive segments
    for _, group_df in cns_df.groupby(is_consecutive.cumsum()):
        # Zip start and end into pairs, then flatten to [start1, end1, start2, end2]
        x_pairs = list(zip(group_df["start"] + f_start_pos(group_df["chrom"]), group_df["end"] + f_start_pos(group_df["chrom"])))
        x = [val for pair in x_pairs for val in pair]
        y = [val for val in group_df[cn_column] for _ in (0, 1)]

        ax.plot(x, y, c=color, linewidth=size, label=label, alpha=alpha, linestyle=linestyle)
        label = None  # only use label for the first segment
    return ax


def plot_dots(ax, cns_df, cn_column, color="green", label=None, alpha=1.0, size=1, assembly=hg19):
    """
    Plots dots representing segments on the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the dots.
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for copy number data.
    color : str, optional
        Color of the dots. Default is "green".
    label : str, optional
        Label for the dots. Default is None.
    alpha : float, optional
        Alpha value for the dots. Default is 1.
    size : float, optional
        Size of the dots. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted dots.
    """
    length = cns_df["end"] - cns_df["start"]
    f_start_pos = _get_start_vector(assembly)
    x = cns_df["start"] + length / 2 + f_start_pos(cns_df["chrom"])
    ax.scatter(x, cns_df[cn_column], s=size, label=label, color=color, alpha=alpha)
    return ax


def plot_bars(ax, cns_df, cn_column, color="green", label=None, alpha=1.0, size=1.0, assembly=hg19):
    """
    Plots bars representing segments on the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the bars.
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for copy number data.
    color : str, optional
        Color of the bars. Default is "green".
    label : str, optional
        Label for the bars. Default is None.
    alpha : float, optional
        Alpha value for the bars. Default is 1.
    size : float, optional
        Line width of the bars. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted bars.
    """
    length = cns_df["end"] - cns_df["start"]
    f_start_pos = _get_start_vector(assembly)
    x = cns_df["start"] + length / 2 + f_start_pos(cns_df["chrom"])
    ax.bar(x, cns_df[cn_column], width=length, color=color, label=label, alpha=alpha, linewidth=size, edgecolor='black')
    return ax


def plot_heatmap(ax, cns_df, cn_column, min_cn = 0, max_cn = 16, assembly=hg19):
    """
    Plots a heatmap of the Copy Number (CN) data.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the heatmap.
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for copy number data.
    min_cn : float, optional
        Minimum copy number value for the color gradient. Default is 0.
    max_cn : int, optional
        Maximum copy number value for the color gradient. Default is 16.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted heatmap.
    """
    f_start_pos = _get_start_vector(assembly)
    # lowest value strictly greater than 0
    f_col = _get_CN_color_vector(min_cn, max_cn)

    ax.set_facecolor("gray")
    labels = []
    for i, ((sample_id), group_df) in enumerate(cns_df.groupby("sample_id")):
        y = -i
        height = 1
        width = group_df["end"] - group_df["start"]
        left = group_df["start"] + f_start_pos(group_df["chrom"])
        color = f_col(group_df[cn_column])
        ax.barh(y, width, height, left=left, color=color)
        labels.append(sample_id)

    ax.set_yticks(-np.arange(len(labels)))
    ax.set_yticklabels(labels)

    ax.autoscale(tight=True)
    return ax


def _get_columns(cns_df, cn_columns):
    if cn_columns == None:
        cn_columns = get_cn_cols(cns_df)
        if len(cn_columns) == 0:
            raise ValueError("If cn_columns is not specified, at least one column matching CN column pattern must be present in data (see documentation for details).")
    elif isinstance(cn_columns, str):
        if not cn_columns in cns_df.columns:
            raise ValueError(f"specified CN column {cn_columns} must in cns_df.columns, have {list(cns_df.columns)}")
        cn_columns = [cn_columns]
    elif isinstance(cn_columns, Sequence):
        if len(cn_columns) <= 0:
            raise ValueError("cn_columns must be a string or a non-empty list of strings")
        elif not all(c in cns_df.columns for c in cn_columns):
            raise ValueError("all elements in cn_columns must be columns in data")
    else:
        raise ValueError("cn_columns must be a string or a list of strings")

    return cn_columns


def _get_colors(colors, line_count):
    if colors == None:
        if line_count == 1:
            colors = ["blue"]
        elif line_count <= 10:
            colors = plt.cm.tab10(np.arange(0, line_count / 10, 1 / 10))
        elif line_count <= 20:
            colors = plt.cm.tab20(np.arange(0, line_count / 20, 1 / 20))
        else:
            colors = plt.cm.hsv(np.linspace(0.0, 1, line_count+1))

    elif line_count == 1:
        colors = [colors]
    elif isinstance(colors, Sequence):
        if line_count != len(colors):
            raise ValueError("colors must be None or a list with the same length as the number of lines")
    else:
        raise ValueError("colors must be None or a list with the same length as the number of lines")
    return colors


def _fig_common(cns_df, f_plot, cn_columns=None, colors=None, size=1, assembly=hg19):
    cn_columns = _get_columns(cns_df, cn_columns)
    groups_df = cns_df.groupby("sample_id")
    line_count = len(groups_df)
    if line_count > 100:
        raise ValueError("Too many samples to plot, please plot fewer than 100 samples at a time (e.g., group by cancer type or other metadata).")
    colors = _get_colors(colors, line_count)
    alpha = (1 / line_count) ** (1 / 3) if f_plot == plot_lines or f_plot == plot_steps else 1 / line_count

    n_columns = len(cn_columns)
    x_min, x_max = x_limits(cns_df, assembly)
    width = max(4, (x_max - x_min) / 200_000_000)
    height = 4 * n_columns
    fig, axes = plt.subplots(n_columns, 1, figsize=(width, height), sharex=True)

    for j, cn_column in enumerate(cn_columns):
        ax = axes[j] if n_columns > 1 else axes

        min_cn = cns_df[cn_column].min()
        max_cn = cns_df[cn_column].max()
        pad = max(0.5, (max_cn - min_cn) * 0.1)
        plot_chr_bg(ax, assembly=assembly, y_min = min_cn - pad, y_max=max_cn + pad, alpha=0.2)
        for i, (group_key, group_df) in enumerate(groups_df):
            color = colors[i]
            label = group_key
            f_plot(
                ax = ax,
                cns_df = group_df,
                cn_column=cn_column,
                color=color,
                label=label,
                alpha=alpha,
                size=size,
                assembly=assembly,
            )

        ax.legend(loc="upper right")
        ax.set_ylabel(f"{cn_column}")
        ax.set_xlim(x_min, x_max)
        if j == n_columns - 1:
            plot_x_ticks(ax, assembly, x_min, x_max)   
            ax.set_xlabel("position on the linear genome")
    
    fig.tight_layout()
    return fig, axes


def fig_lines(cns_df, cn_columns=None, colors=None, size=1, assembly=hg19):
    """
    Creates a line plot for each of the CN columns.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    colors : list of str, optional
        List of colors to use for the plots. If None, colors are generated automatically.
    size : int, optional
        Size of the feature of the plot - line/boundary width or dot size. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    list of matplotlib.axes.Axes
        List of axes in the figure.
    """
    return _fig_common(cns_df, plot_lines, cn_columns, colors, size, assembly)


def fig_dots(cns_df, cn_columns=None, colors=None, size=1, assembly=hg19):
    """
    Creates a dot plot for each of the CN columns.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    colors : list of str, optional
        List of colors to use for the plots. If None, colors are generated automatically.
    size : int, optional
        Size of the plot. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    list of matplotlib.axes.Axes
        List of axes in the figure.
    """
    return _fig_common(cns_df, plot_dots, cn_columns, colors, size, assembly)


def fig_bars(cns_df, cn_columns=None, colors=None, size=1, assembly=hg19):
    """
    Creates a bar plot for each of the CN columns.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    colors : list of str, optional
        List of colors to use for the plots. If None, colors are generated automatically.
    size : int, optional
        Size of the plot. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    list of matplotlib.axes.Axes
        List of axes in the figure.
    """
    return _fig_common(cns_df, plot_bars, cn_columns, colors, size, assembly)
 

def fig_steps(cns_df, cn_columns=None, colors=None, size=1, assembly=hg19):
    """
    Creates a step plot for each of the CN columns.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    colors : list of str, optional
        List of colors to use for the plots. If None, colors are generated automatically.
    size : int, optional
        Size of the plot. Default is 1.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    list of matplotlib.axes.Axes
        List of axes in the figure.
    """
    return _fig_common(cns_df, plot_steps, cn_columns, colors, size, assembly)
 


def _make_layout(width, height, n_columns, vertical):
    if vertical:
        height = height * n_columns
        sharex, sharey = True, False
        n_rows, n_cols = n_columns, 1
    else:
        width = width * n_columns
        sharex, sharey = False, True
        n_rows, n_cols = 1, n_columns
    width = max(3, width)
    height = max(3, height)
    return plt.subplots(n_rows, n_cols, figsize=(width, height), sharex=sharex, sharey=sharey)


def fig_heatmap(cns_df, cn_columns=None, min_cn = 0, max_cn = 10, vertical = None, assembly=hg19):
    """
    Creates a heatmap figure from copy number segment (CNS) data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data with columns for sample_id, chromosome positions,
        and copy number values.
    cn_columns : list of str, optional
        List of column names containing copy number values. If None, columns are
        inferred from the DataFrame.
    min_cn : int, optional
        Minimum copy number value for the color scale. Default is 0.
    max_cn : int, optional
        Maximum copy number value for the color scale. Default is 10.
    vertical : bool, optional
        Whether to stack plots vertically. If None, orientation is determined by
        the aspect ratio. Default is None.
    assembly : Assembly object, optional
        Genome assembly object defining chromosome sizes. Default is hg19.

    Returns
    -------
    tuple
        A tuple containing:
        - matplotlib.figure.Figure: The created figure
        - numpy.ndarray: Array of matplotlib.axes.Axes objects
    
    Notes
    -----
    The figure size is automatically determined based on the genome size and 
    number of samples. The layout (vertical vs horizontal) is determined by
    the aspect ratio unless explicitly specified.
    """
    cn_columns = _get_columns(cns_df, cn_columns)

    sample_count = len(cns_df["sample_id"].unique())
    n_columns = len(cn_columns)
    x_min, x_max = x_limits(cns_df, assembly)
    width = (x_max - x_min) / 200_000_000
    height = sample_count / 5
    vertical = vertical if vertical != None else width > height
    fig, axes = _make_layout(width, height, n_columns, vertical)
    
    min_cn = max(cns_df[cn_columns][cns_df[cn_columns] > 0].min().min(), min_cn)
    max_cn = min(cns_df[cn_columns].max().max(), max_cn)

    for j, cn_column in enumerate(cn_columns):
        ax = axes[j] if n_columns > 1 else axes
        plot_heatmap(
            ax = ax,
            cns_df = cns_df,
            cn_column=cn_column,
            min_cn = min_cn,
            max_cn = max_cn,
            assembly=assembly,
        )
        if vertical:
            ax.set_ylabel(f"{cn_column}")
        elif j == 0:
            ax.set_ylabel("")
        if not vertical:
            ax.set_xlabel(f"{cn_column}")
        elif j == n_columns - 1:
            ax.set_xlabel(f"position on the linear genome")
        plot_x_lines(ax, assembly)
        plot_x_ticks(ax, assembly, x_min, x_max)   
        ax.set_xlim(x_min, x_max)
    
    ax.margins(x=0, y=0)

    # Add legend
    handles = []
    handles.append(mpatches.Patch(facecolor ='blue', label=f'{max_cn:.2f}', edgecolor='black'))
    handles.append(mpatches.Patch(facecolor ='white', label=f'{min_cn:.2f}', edgecolor='black'))
    handles.append(mpatches.Patch(facecolor ='red', label='0', edgecolor='black'))
    handles.append(mpatches.Patch(facecolor ='gray', label='NaN', edgecolor='black'))
    last_ax = axes[-1] if n_columns > 1 else axes
    last_ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1))
    
    return fig, axes


def x_limits(cns_df, assembly=hg19):
    """
    Calculates the x-axis limits for the plot based on the CNS data and genome assembly.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    tuple
        A tuple containing the minimum and maximum x-axis limits.
    """
    offset = _get_start_vector(assembly)(cns_df["chrom"])
    min_x = (cns_df["start"] + offset).min()
    max_x = (cns_df["end"] + offset).max()
    return min_x, max_x


def y_limits(cns_df, column):
    """
    Get the limits for the y-axis based on the CNS data and the column.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    column : str
        Column name for the y-axis data.

    Returns
    -------
    tuple
        A tuple containing the minimum and maximum y-axis limits.

    """
    min_y = cns_df[column].min()
    max_y = cns_df[column].max()
    return min_y, max_y


def _plot_rectangles(ax, items, y_min, y_max, assembly, color_func, alpha):
    height = y_max - y_min

    for item in items:
        chrom, start, end = item[:3]
        x_pos = start + assembly.chr_starts[chrom]
        color = color_func(item)
        width = end - start
        rect = mpatches.Rectangle((x_pos, y_min), width, height, color=color, alpha=alpha)
        ax.add_patch(rect)

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, assembly.gen_len)


def plot_chr_bg(ax, y_min=0, y_max=2, assembly=hg19, alpha=0.2):
    """
    Plots the chromosome background on the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the chromosome background.
    assembly : object, optional
        Genome assembly to use. Default is hg19.
    y_min : float, optional
        Minimum y-axis value. Default is 0.
    y_max : float, optional
        Maximum y-axis value. Default is 1.
    alpha : float, optional
        Alpha value for the background. Default is 0.2.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted chromosome background.
    """
    def color_func(item):
        chrom = item[0]
        is_even = item[1] % 2 == 0
        return "darkgray" if is_even else "lightgray"

    items = [(chrom, i, i + length) for i, (chrom, length) in enumerate(assembly.chr_lens.items())]
    _plot_rectangles(ax, items, y_min, y_max, assembly, color_func, alpha=alpha)


def plot_cytobands(ax, y_min=0, y_max=2, assembly=hg19, alpha=0.2, color=None):
    """
    Plots cytobands on the background of the ax.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the cytobands.
    bounds : tuple, optional
        Bounds for the plot (x_min, y_min, x_max, y_max). If None, it is inferred from the data.
    assembly : object, optional
        Genome assembly to use. Default is hg19.
    alpha : float, optional
        Alpha value for the cytobands. Default is 0.2.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted cytobands.
    """
    f_color = lambda item: color if color is not None else cytoband_color[item[4]]
    _plot_rectangles(ax, assembly.cytobands, y_min, y_max, assembly, f_color, alpha)


def plot_gaps(ax, y_min=0, y_max=2, assembly=hg19, alpha=0.2, color=None):
    f_color = lambda item: color if color is not None else gap_color[item[3]]
    _plot_rectangles(ax, assembly.gaps, y_min, y_max, assembly, f_color, alpha)


def _create_label(num):
    if num < 1_000:
        return str(num)
    elif num < 1_000_000:
        return f"{num // 1_000}kb"
    else:
        return f"{num // 1_000_000}mb"

def plot_x_ticks(ax, assembly=hg19, min_x=0, max_x=None):
    """
    Plots the x-axis ticks for the chromosomes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the x-axis ticks.
    assembly : object, optional
        Genome assembly to use. Default is hg19.
    x_min : float, optional
        Minimum x-axis value. If None, it is inferred from the data.
    x_max : float, optional
        Maximum x-axis value. If None, it is inferred from the data.

    Returns
    -------
    list of float
        List of x-axis tick positions.
    """
    positions = list(assembly.chr_lens.items())
    if max_x is None:
        max_x = assembly.gen_len

    x_pos = 0
    major_tick_pos = []
    minor_tick_pos = []
    minor_tick_labels = []
    for chrom, length in positions:
        if min_x <= x_pos <= max_x:
            label_text = "\n" + chrom[3:]
            major_tick_pos.append(x_pos)
            minor_tick_pos.append(x_pos + length / 2)
            minor_tick_labels.append(label_text)
        # Update the x position for the next chromosome
        x_pos += length

    if len(major_tick_pos) < 2:
        ax.set_xticks([min_x, max_x])
        ax.set_xticklabels([_create_label(min_x), _create_label(max_x)], rotation=90)
        return
    ax.set_xticks(major_tick_pos)
    ax.set_xticks(minor_tick_pos, minor=True)
    ax.set_xticklabels([" "] * len(major_tick_pos))
    ax.set_xticklabels(minor_tick_labels, minor=True)
    ax.set_xlim(min_x, max_x)
        # Hide the lines for the minor ticks
    ax.tick_params(axis='x', which='minor', length=0, pad=-8, labelsize=10)

    return major_tick_pos, minor_tick_pos


def add_cytoband_legend(ax):
    """
    Adds a legend for cytobands to the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to which the legend is added.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the added legend.
    """
    legend_elements = [mpatches.Patch(color=color, label=gie_stain) for gie_stain, color in cytoband_color.items()]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc="upper left")


def add_gap_legend(ax):
    """
    Adds a legend for gaps to the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to which the legend is added.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the added legend.
    """
    legend_elements = [mpatches.Patch(color=color, label=gap_type) for gap_type, color in gap_color.items()]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc="upper left")


def no_y_ticks(ax):
    """
    Removes y-axis ticks from the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes from which to remove the y-axis ticks.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the y-axis ticks removed.
    """
    ax.set_yticks([])
    ax.set_yticklabels([])


def no_x_ticks(ax):
    """
    Removes x-axis ticks from the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes from which to remove the x-axis ticks.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the x-axis ticks removed.
    """
    ax.set_xticks([])
    ax.set_xticklabels([])


def plot_x_lines(ax, assembly=hg19, positions=None, width=1, alpha=0.5):
    """
    Plots vertical lines at chromosome boundaries on the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes on which to plot the vertical lines.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted vertical lines.
    """
    positions = list(accumulate(assembly.chr_lens.values()))
    for pos in positions:
        ax.axvline(
            pos,
            color="black",
            linewidth=width,
            linestyle="--",
            alpha=alpha,
        )
