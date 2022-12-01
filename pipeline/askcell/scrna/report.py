"""Report result for generalized scRNA-seq analysis."""

from __future__ import annotations

from datetime import date

from fpdf import FPDF

from .. import Path, PathLike
from ..utils.types import ImmutableStrMapping


def write_report(
    outputs: ImmutableStrMapping,
    params: ImmutableStrMapping,
    report: PathLike,
    runid: str,
    orientation: str = "L",
    unit: str = "mm",
    format: str = "A4",
    font: str = "Arial",
    size: int = 16,
    style: str = "B",
) -> None:
    """Write report.

    Args:
        outputs: mapping of key word and output filenames.
        report: report filename

    """
    # 1. Set up the PDF doc basics
    pdf = FPDF(orientation=orientation, unit=unit, format=format)
    pdf.add_page()
    pdf.set_font(font, style=style, size=size)

    # 2. Layout the PDF doc contents
    pdf.cell(40, 10, f'{runid} report, {date.today().strftime("%b-%d-%Y")}')
    pdf.ln(10)

    pdf.set_font("Arial", style="", size=12)
    pdf.cell(txt="Following parameters are used during the analysis:")
    pdf.ln(5)
    pdf.set_font("Arial", style="", size=8)
    for k, v in params.items():
        if not isinstance(v, (Path, dict)):
            pdf.cell(txt=f"* {k} = {v}")
            pdf.ln(5)

    for k, v in params[params["chemistry"]].items():
        if not isinstance(v, (Path, dict)):
            pdf.cell(txt=f"* {k} = {v}")
            pdf.ln(5)

    pdf.ln(5)
    pdf.cell(txt="Following parameters are used during PCA analysis:")
    pdf.ln(5)
    for k, v in params["pca"].items():
        if not isinstance(v, (Path, dict)):
            pdf.cell(txt=f"* {k} = {v}")
            pdf.ln(5)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.ln(10)
    pdf.cell(txt="Violin Plots")
    pdf.ln(10)
    pdf.set_font("Arial", style="", size=12)
    pdf.cell(txt=", ".join(params["preQC"]["violin_plots"]["keys"]))
    pdf.ln(5)
    pdf.cell(0, 10, "Violin Plot - pre-filtering", 0, 1)
    pdf.image(outputs["violin_plot_pre"], w=pdf.epw)

    pdf.ln(5)
    pdf.cell(0, 10, "Violin Plot - post-filtering", 0, 1)
    pdf.image(outputs["violin_plot"], w=pdf.epw)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Scatter Plots")
    pdf.ln(5)
    pdf.set_font("Arial", style="", size=12)
    pdf.cell(0, 10, "Scatter Plot - pre- and post-filtering", 0, 1)
    pdf.image(outputs["scatter_plot_pre"], w=pdf.epw / 2, y=pdf.eph * 0.1)
    pdf.set_y(0)
    pdf.image(outputs["scatter_plot"], w=pdf.epw / 2, x=pdf.epw / 2, y=pdf.eph * 0.1)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Box Plots of highly expressed genes - pre- and post-filtering")
    pdf.ln(5)
    pdf.set_font("Arial", style="", size=12)
    pdf.cell(0, 10, "Highly expressed genes - pre- and post-filtering", 0, 1)
    pdf.image(outputs["box_plot_pre"], w=pdf.epw / 2, y=pdf.eph * 0.1)
    pdf.set_y(0)
    pdf.image(outputs["box_plot"], w=pdf.epw / 2, x=pdf.epw / 2, y=pdf.eph * 0.1)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Sex-bias Plots")
    pdf.ln(5)
    pdf.set_font("Arial", style="", size=12)
    pdf.cell(0, 10, "Presence of chrY expressions and XIST-genes", 0, 1)
    pdf.image(outputs["sexbias_plot1"], w=pdf.epw / 2, y=pdf.eph * 0.1)
    pdf.set_y(0)
    pdf.image(outputs["sexbias_plot2"], w=pdf.epw / 2, x=pdf.epw / 2, y=pdf.eph * 0.1)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Additional Plots")
    pdf.ln(5)
    pdf.set_font("Arial", style="", size=12)
    pdf.cell(0, 10, "S phase and G2M phase scores", 0, 1)
    pdf.image(outputs["cellcycle_plot"], w=pdf.epw * 0.9)
    pdf.ln(5)
    pdf.cell(0, 10, "Doublet prediction plots", 0, 1)
    pdf.image(outputs["doublet_plot"], w=pdf.epw * 0.8)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Dimensional Reduction Plots - Highly Variable Genes")
    pdf.ln(5)
    pdf.image(outputs["hvg_plot"], w=pdf.epw * 0.9)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Dimensional Reduction Plots - PCs")
    pdf.ln(5)
    pdf.image(outputs["pca_plot"], h=pdf.eph * 0.9)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Dimensional Reduction Plots - UMAP")
    pdf.ln(5)
    pdf.image(outputs["umap_plot"], w=pdf.epw * 0.90)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Dimensional Reduction Plots - Rank Genes Groups")
    pdf.ln(5)
    pdf.image(outputs["rank_plot"], h=pdf.eph * 0.95)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Marker Expression Plots")
    pdf.ln(5)
    pdf.image(outputs["marker_umap_plot"], h=pdf.eph * 0.95)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Marker Expression Plots")
    pdf.ln(5)
    pdf.image(outputs["marker_dotplot"], h=pdf.eph * 0.95)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Marker Expression Plots")
    pdf.ln(5)
    pdf.image(outputs["marker_tracks_plot"], w=pdf.epw * 0.9, h=pdf.eph * 0.95)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Marker Expression Plots")
    pdf.ln(5)
    pdf.image(outputs["marker_matrix_plot"], h=pdf.eph * 0.95)

    pdf.add_page()
    pdf.set_font("Arial", style="B", size=16)
    pdf.cell(txt="Cluster Proportion Plots")
    pdf.ln(5)
    pdf.image(outputs["props_plot"], h=pdf.eph * 0.95)

    pdf.output(report, "F")
