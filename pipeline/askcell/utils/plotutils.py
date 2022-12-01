"""plots.py for various plots."""

from __future__ import annotations

import pandas as pd
import plotnine as p9


def facet(
    row: str | None,
    col: str | None,
    nrow: int | None = 1,
    ncol: int | None = 1,
) -> p9.facets.facet_wrap.facet_wrap | p9.facets.facet_grid.facet_grid:
    """Facet aggregated plot layout for creating multi-panel plots.

    Args:
        row: row variable to create multi-panel plots, facet_grid is used if row and col are used
        col: column variable to create multi-panel plots, facet_grid is used if row and col are used
        nrow: number of rows
        ncol: number of columns

    Returns:
        facet_grid or facet_wrap for multi-panel plots or default theme

    """
    if row and col:
        return p9.facet_grid(f"{row} ~ {col}")
    elif row:
        return p9.facet_wrap(row, nrow=nrow)
    elif col:
        return p9.facet_wrap(col, ncol=ncol)
    else:
        return p9.theme()


def plot_histogram(
    data: pd.DataFrame,
    *,
    x: str,
    binwidth: int = 25,
    row: str | None = None,
    col: str | None = None,
    nrow: int | None = None,
    ncol: int | None = None,
    fill: str | None = None,
    title: str = "",
    xlabel: str = "",
    xlabel_rotation: int = 90,
    axis_text_size: int | None = None,
    width: int | None = 4,
    height: int | None = 3,
    dpi: int = 100,
) -> p9.ggplot.ggplot:
    """Plot histogram.

    Args:
        data: dataframe containing data to plot
        x: variable to be used for x-axis
        binwidth: bin width for histogram
        row: row variable used for facet_wrap, facet_grid is used if row and col are used
        col: column variable used for facet_wrap, facet_grid is used if row and col are used
        nrow: number of rows used for facetting
        ncol: number of columns used for facetting
        fill: if True, fill in the box by variable values
        title: title of the plot
        x_label: set x-axis label
        xlabel_rotation: set the rotation of the text
        axis_text_size: customize the size of the text on the axis labelling
        width: Figure width in inches
        height: Figure height in inches
        dpi: Figure DPI resolution

    Returns:
        plotnine histogram object

    """
    return (
        p9.ggplot(data, p9.aes(x=x, fill=fill or x))
        + p9.geom_histogram(binwidth=binwidth)
        + p9.ggtitle(title)
        + p9.xlab(xlabel)
        + p9.theme(
            axis_text_x=p9.element_text(rotation=xlabel_rotation, vjust=1, hjust=1, size=axis_text_size),
            figure_size=(width, height),
            dpi=dpi,
        )
        + facet(row=row, col=col, nrow=nrow, ncol=ncol)
    )


def plot_boxplot(
    data: pd.DataFrame,
    *,
    x: str,
    y: str,
    row: str | None = None,
    col: str | None = None,
    nrow: int | None = None,
    ncol: int | None = None,
    fill: str | None = None,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    xlabel_rotation: int = 90,
    axis_text_size: int | None = None,
    show_legend: bool = False,
    legend_position: str = "top",
    points_position: str = "jitterdodge",
    width: int | None = 4,
    height: int | None = 3,
    dpi: int = 100,
) -> p9.ggplot.ggplot:
    """Boxplot counts.

    Args:
        data: dataframe containing data to plot
        x: variable to be used for x-axis
        y: variable to be used for y-axis
        row: row variable used for facet_wrap, facet_grid is used if row and col are used
        col: column variable used for facet_wrap, facet grid is used if both row and col are used
        fill: if True, fill in the box by variable values
        title: title of the plot
        x_label: set x-axis label
        y_label: set y-axis label
        xlabel_rotation: set the rotation of the text
        axis_text_size: customize the size of the text on the axis labelling
        show_legend: flag to show legend, default is False
        legend_position: legend position
        points_position: boxplot point position: identity or jitterdodge, default is jitterdodge
        width: Figure width in inches
        height: Figure height in inches
        dpi: Figure DPI resolution

    Returns:
        plotnine boxplot object

    """
    return (
        p9.ggplot(data, p9.aes(x=x, y=y, fill=fill or x))
        + p9.geom_boxplot(show_legend=show_legend)
        + p9.ggtitle(title)
        + p9.xlab(xlabel)
        + p9.ylab(ylabel)
        + p9.theme(
            axis_text_x=p9.element_text(rotation=xlabel_rotation, vjust=1, hjust=1, size=axis_text_size),
            legend_position=legend_position,
            figure_size=(width, height),
            dpi=dpi,
        )
        + facet(row=row, col=col, nrow=nrow, ncol=ncol)
    )


def plot_points(
    data: pd.DataFrame,
    *,
    x: str,
    y: str,
    row: str | None = None,
    col: str | None = None,
    nrow: int | None = None,
    ncol: int | None = None,
    color: str | None = None,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    xlabel_rotation: int = 90,
    axis_text_size: int | None = None,
    size: int = 2,
    alpha: float | None = 0.5,
    width: int | None = 4,
    height: int | None = 3,
    position: str | None = "identity",
    blank_x: bool = False,
    dpi: int = 100,
) -> p9.ggplot.ggplot:
    """Plot data points.

    Args:
        data: dataframe containing data to plot
        x: variable to be used for x-axis
        y: variable to be used for y-axis
        row: row variable used for facet_wrap, facet_grid is used if row and col are used
        col: column variable used for facet_wrap, facet grid is used if both row and col are used
        nrow: number of rows used for facetting
        ncol: number of columns used for facetting
        color: color of the points on the plot
        title: title of the plot
        x_label: set x-axis label
        y_label: set y-axis label
        xlabel_rotation: set the rotation of the text
        axis_text_size: customize the size of the text on the axis labelling
        size: size of the point on the plot
        alpha: transparency of the point on the plot
        width: Figure width in inches
        height: Figure height in inches
        position: position adjustment
        blank_x: if set True, leave x-axis tick and text blank
        dpi: Figure DPI resolution

    Returns:
        plotnine plot object

    """
    return (
        p9.ggplot(data, p9.aes(x=x, y=y, color=color))
        + p9.geom_point(size=size, alpha=alpha, position=position)
        + p9.ggtitle(title)
        + p9.xlab(xlabel)
        + p9.ylab(ylabel)
        + p9.theme(
            figure_size=(width, height),
            dpi=dpi,
            axis_ticks_major_x=p9.element_blank() if blank_x else None,
            axis_text_x=p9.element_blank()
            if blank_x
            else p9.element_text(
                rotation=xlabel_rotation,
                vjust=1,
                hjust=1,
                size=axis_text_size,
            ),
        )
        + facet(row=row, col=col, nrow=nrow, ncol=ncol)
    )
