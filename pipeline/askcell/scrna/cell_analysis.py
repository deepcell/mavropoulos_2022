"""Cell-level analysis task for scRNA-seq.

Perform downstream analysis using Scanpy.
This is currently focused around 10x platform data analysis based on
`Current best practices in single-cell RNA-seq analysis: a tutorial`
by Malte D Luecken and Fabian J thesis (2019).

After aggregating the raw counts from multiple samples,
it applies the normalization and transformation, along with QC analysis
and provides dimensionally reduced data in AnnData or Seurat format.

These are the steps in this task.

    * Aggregate counts into a single data
    * Performs multiple operation to the data using scanpy
        - QC and filtering of poor quality cells and genes
        - normalization / transformation
        - PCa
        - association with the marker genes
    * Collect and write metrics
    * Copy local results files to final output destination
    * Remove intermediate files


Example command:
    scrna --local-temp-dir /home/julie/tmp/ analyze \
        --input-dir /home/julie/tmp/test_scrna \
        --output-dir /home/julie/tmp/test_scrna/scanpy \
        --threads 8 /home/julie/tools/pipeline/tests/data/test_scrna/request.xlsx \
        /home/julie/tools/pipeline/tests/data/test_scrna/samples2.tsv

"""

from __future__ import annotations

import logging

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import seaborn as sns

from matplotlib import pyplot as plt

from .. import PathLike
from ..utils.config import PackageConfig
from ..utils.types import ImmutableStrMapping
from .report import write_report
from .request import AnalysisRequest


def run_scanpy_analysis(
    inputs: ImmutableStrMapping,
    outputs: ImmutableStrMapping,
    scratch: ImmutableStrMapping,
    params: ImmutableStrMapping,
    request: AnalysisRequest,
    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]],
    **kwargs: Any,
) -> None:
    """Run scanpy-based downstream analysis as a step."""
    adata = load_counts({sid: sample.data for sid, (_, sample) in sample_configs.items()})
    # FIXME: Make figure settings configurabe
    sc.settings.set_figure_params(dpi=200, figsize=[12, 8], fontsize=15)

    logging.info("[Task] Performing pre-QC.")
    perform_pre_qc(adata, params).to_csv(scratch["preqc_summary"], sep="\t", index=False)
    plot_qc(adata, params["preQC"])

    logging.info("[Task] Filtering by min genes and cells.")
    sc.pp.filter_cells(adata, min_genes=params["min_genes"])
    sc.pp.filter_genes(adata, min_cells=params["min_cells"])
    logging.info(f" ====> Remaining cells {adata.n_obs}, genes {adata.n_vars}")

    logging.info("[Task] Save non-sliced aggregated data.")
    adata.write_h5ad(scratch["raw_aggregated_h5ad"])

    logging.info("[Task] Filtering cells by QC, this task slices the data.")
    adata = slice_filtered_cells(adata, params)
    logging.info(f" ====> Remaining cells {adata.n_obs}, genes {adata.n_vars}")

    logging.info("[Task] Filtering genes by QC, this task slices the data.")
    adata = slice_filtered_genes(adata, params)
    logging.info(f" ====> Remaining cells {adata.n_obs}, genes {adata.n_vars}")

    logging.info("[Task] QC sex bias.")
    qc_sex_bias(adata, params)

    logging.info("[Task] Plot QC-checked data.")
    plot_qc(adata, params["QC"])
    plot_sexbias(adata, params["QC"])

    logging.info("[Task] Check doublets.")
    # revert back to the raw counts as the main matrix in adata
    adata = check_doublets(adata, params["doublets"])

    logging.info("[Task] Normalize data.")
    normalize(
        adata,
        target_depth=params["normalize"]["target_depth"],
        is_to_log=params["normalize"]["is_to_log"],
        is_to_scale=params["normalize"]["is_to_scale"],
    )

    # set .raw attribute to the normalized for later use
    logging.info("[Task] Save the normalized counts in the raw slot.")
    adata.raw = adata

    logging.info("[Task] Save filtered and sliced raw and normalized data.")
    adata.write_h5ad(scratch["qc_filtered_h5ad"])

    logging.info("[Task] Check cell cycle genes.")
    score_cellcycle_phase(
        adata,
        [x.strip() for x in open(inputs["cellcycle_genes"])],
        params,
    )

    sc.settings.set_figure_params(dpi=200, figsize=[12, 8], fontsize=10)
    sc.pl.violin(
        adata,
        params["QC"]["cellcycle_plots"]["keys"],
        jitter=0.4,
        groupby=params["QC"]["cellcycle_plots"]["groupby"],
        rotation=30,
        show=False,
        use_raw=False,
    )
    plt.savefig(params["QC"]["cellcycle_plots"]["filename"])
    plt.close()

    if params["correct_batch"]:
        adata = correct_batch(adata, key="runid")

    logging.info("[Task] Perform dimensional reduction.")
    adata = reduce_dimensionality(adata, **params["pca"])

    logging.info("[Task] Save a table with top ranked genes in a group.")
    get_top_ranked_genes(adata, **params["gene_ranks"])

    logging.info("[Task] Get cluster proportions and plot.")
    props_param = params["cluster_proportion"]
    plot_cluster_proportions(
        get_cluster_proportions(
            adata,
            cluster_key=props_param["cluster_key"],
            sample_key=props_param["sample_key"],
        ),
        props_param["filename"],
    )

    logging.info("[Task] Show clusters expressing the markers of interest.")
    plot_marker_expressions(adata, params["marker_plots"])

    logging.info("[Task] Save normalized data with PCA embeddings.")
    adata.write_h5ad(scratch["h5ad"])

    write_report(scratch, params, scratch["pdf_report"], runid=params["requestid"])


def load_counts(raw_counts: dict[str, ImmutableStrMapping]) -> ad.AnnData:
    """Load count metrices.

    Args:
        raw_counts: key and value pair of sid and raw count matrices paths, label

    Returns:
        counts data in AnnData object

    """
    data = []
    for sid, sample in raw_counts.items():
        d = sc.read_10x_h5(sample["input"]["h5_count"])
        d.var_names_make_unique()
        d.obs["sample"] = sid
        d.obs["type"] = sample["label"] or "sample"
        d.obs["runid"] = sample["runid"]
        data.append(d)

    return ad.AnnData.concatenate(*data)


def perform_pre_qc(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> pd.DataFrame:
    """Perform preQC on count data.

    Calculate the percentage of mitocondrial and ribosomal genes per cell.

    Citing from "Simple Single Cell" workflows (Lun, McCarthy & Marioni, 2017):
    "High proportions are indicative of poor-quality cells
    (Islam et al. 2014; Ilicic et al. 2016),
    possibly because of loss of cytoplasmic RNA from perforated cells.
    The reasoning is that mitochondria are larger than individual transcript molecules
    and less likely to escape through tears in the cell membrane."

    Args:
        adata: aggregated raw counts data
        params: QC parameters

    Returns:
        metrics data on sample count performance

    """
    # Define which genes are mitochondrial, ribosomal and hemoglogin
    adata.var["mt"] = adata.var_names.str.startswith(params["mt_gene_name"])
    adata.var["ribo"] = adata.var_names.str.startswith(params["ribo_gene_name"])
    adata.var["hb"] = adata.var_names.str.contains(params["hb_gene_name_pattern"])

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    mito_genes = adata.var_names.str.startswith(params["mt_gene_name"])
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs["percent_mt2"] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1

    summary = adata.obs.groupby("sample").describe()
    summary.columns = ["_".join(col) for col in summary.columns.values]
    return summary.reset_index()


def collect_cells_metrics() -> None:
    """Collect sample-level cells metrics.

    Args:
        result_dir: path to alignment/counts data

    Returns:
        metrics data on sample alignment and count performance

    """
    pass


def slice_filtered_cells(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> ad.AnnData:
    """Filter out cells by applying QC, including doublets.

    Args:
        adata: aggregated raw counts data
        params: thresholds needed for filtering step

    Returns:
        sliced Ann Data

    """
    chemistry = params["chemistry"]

    keep = (
        (adata.obs["n_genes_by_counts"] < params[chemistry]["max_n_genes_by_counts"])
        & (adata.obs["n_genes_by_counts"] > params[chemistry]["min_n_genes_by_counts"])
        & (adata.obs["pct_counts_mt"] < params["max_pct_counts_mt"])
        & (adata.obs["pct_counts_ribo"] > params["min_pct_counts_ribo"])
        & (adata.obs["pct_counts_ribo"] < params["max_pct_counts_ribo"])
    )

    return adata[keep, :]


def slice_filtered_genes(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> ad.AnnData:
    """Filter out genes by applying QC.

    Args:
        adata: aggregated raw counts data
        params: thresholds needed for filtering step

    Returns:
        sliced Ann Data

    """
    if not adata.var_names.size:
        # no data left to filter
        return adata

    remove = np.array([False] * adata.var_names.size)

    for gene in params["filter_genes"]:
        remove = np.add(remove, adata.var_names.str.startswith(gene))

    if params["remove_mito_genes"]:
        remove = np.add(remove, adata.var_names.str.startswith(params["mt_gene_name"]))

    if params["remove_hb_genes"]:
        remove = np.add(remove, adata.var_names.str.contains(params["hb_gene_name_pattern"]))

    return adata[:, np.invert(remove)]


def qc_sex_bias(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> None:
    """QC sex bias.

    Identify reads from chrY (males) and XIST (mainly females) to determine the sex.
    And detect any sample mixups, if the sample metadata sex does not agree.

    Args:
        adata: aggregated raw counts data
        params: thresholds needed for sex bias correction

    """
    annot = sc.queries.biomart_annotations(
        "hsapiens",
        [
            "ensembl_gene_id",
            "external_gene_name",
            "start_position",
            "end_position",
            "chromosome_name",
        ],
    ).set_index("external_gene_name")

    chrY_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
    adata.obs["percent_chrY"] = np.sum(adata[:, chrY_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100

    if adata.X[:, adata.var_names.str.match("XIST")].size:
        adata.obs["XIST-counts"] = adata.X[:, adata.var_names.str.match("XIST")].toarray()
    else:
        adata.obs["XIST-counts"] = np.zeros(adata.obs_names.size)


def score_cellcycle_phase(
    adata: ad.AnnData,
    genes: list[str],
    params: ImmutableStrMapping,
) -> None:
    """Score cellcycle phase.

    Update data with a score for S phase and G2M phase and the predicted cell cycle phase.

    Args:
        adata: aggregated raw counts data
        genes: cell cycle genes
        params: parameters for cellcycle score calculation

    """
    # FIXME: remove 43! update cellcycle metadata with phase information column
    s_genes = genes[:43]
    g2m_genes = genes[43:]

    # cellcycle_genes = [x for x in genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)


def normalize(
    adata: ad.AnnData,
    target_depth: int = 10000,
    is_to_log: bool = True,
    is_to_scale: bool = True,
) -> None:
    """Normalize and transform data.

    Args:
        adata: aggregated raw counts data
        is_to_log: if set true, logaritmize
        is_to_scale: if set true, transform

    """
    # normalize to depth 10 000
    sc.pp.normalize_total(adata, target_sum=target_depth)

    # logaritmize
    if is_to_log:
        sc.pp.log1p(adata)

    if is_to_scale:
        sc.pp.scale(adata)


def reduce_dimensionality(
    adata: ad.AnnData,
    *,
    pca_plot: PathLike,
    umap_plot: PathLike,
    hvg_plot: PathLike,
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
    max_sd: float = 10,
    n_pcs: int = 30,
    n_neighbors: int = 20,
    n_top_genes: int = 5000,
    svd_solver: str = "arpack",
    use_highly_variable_genes_only: bool = False,
    regress_cellcycle_genes: bool = False,
    regress_out_variables: list[str] | None = None,
) -> ad.AnnData:
    """Perform dimensional reduction.

    Args:
        adata: aggregated counts data
        filename: umap plot filename
        FIXME: add additional parameters needed for dimensional reduction

    Returns:
        Ann Data with dimensional reduction

    """
    # compute variable genes
    sc.pp.highly_variable_genes(
        adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        batch_key="runid",
        n_top_genes=n_top_genes,
    )
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(hvg_plot)
    plt.close()
    logging.info(f" ====> Highly variable genes: {sum(adata.var.highly_variable)}")

    # subset for variable genes in the dataset, not needed for PCA as it auto-detects
    # if use_highly_variable_genes_only:
    #    adata = adata[:, adata.var["highly_variable"]]
    # adata = adata.copy() #run this line if you get the "AttributeError: swapaxes not found"

    # regress out unwanted variables
    if regress_cellcycle_genes and regress_out_variables is not None:
        regress_out_variables += ["S_score", "G2M_score"]
    if regress_out_variables is not None:
        sc.pp.regress_out(adata, keys=list(regress_out_variables))

    # scale data, clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=max_sd)
    sc.tl.pca(adata, svd_solver=svd_solver, n_comps=n_pcs, use_highly_variable=use_highly_variable_genes_only)

    # FIXME: separate out the plotting function
    sc.settings.set_figure_params(dpi=200, figsize=[8, 8], fontsize=15)
    sc.pl.pca_overview(
        adata,
        color=["sample", "type"],
        components=["1,2", "2,3"],
        show=False,
    )
    plt.savefig(pca_plot)
    plt.close()

    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.umap(adata)

    # run leiden clustering, it directly clusters neighbors graph of cells
    sc.tl.leiden(adata, key_added="leiden_1.0")  # default resolution in 1.0
    sc.tl.leiden(adata, resolution=0.6, key_added="leiden_0.6")
    sc.tl.leiden(adata, resolution=0.4, key_added="leiden_0.4")
    sc.tl.leiden(adata, resolution=1.4, key_added="leiden_1.4")

    # FIXME: separate out the plotting function
    sc.pl.umap(
        adata,
        color=[
            "leiden_1.0",
            "leiden_0.6",
            "leiden_0.4",
            "leiden_1.4",
            "sample",
            "type",
            "doublet_scores",
            "S_score",
            "G2M_score",
        ],
        show=False,
    )
    plt.savefig(umap_plot)
    plt.close()

    return adata


def get_top_ranked_genes(
    adata: ad.AnnData,
    *,
    rank_group: str,
    num_genes_from_rank: int,
    filename: PathLike,
    rank_plot: PathLike,
    method: str = "wilcoxon",
) -> None:
    """Get top ranked genes.

    Args:
        adata: aggregated counts data
        filename: filename to save ranked genes
        FIXME: add additional parameters needed for dimensional reduction

    Returns:
        Ann Data with dimensional reduction

    """
    sc.tl.rank_genes_groups(adata, rank_group, method=method)
    exp_result = adata.uns["rank_genes_groups"]

    pd.DataFrame(
        {
            group + "_" + key: exp_result[key][group]
            for group in exp_result["names"].dtype.names
            for key in ["names", "pvals"]
        },
    ).to_csv(filename, sep="\t", index=True)

    sc.pl.rank_genes_groups(adata, n_genes=num_genes_from_rank, show=False, sharey=False)
    plt.savefig(rank_plot)
    plt.close()


def plot_marker_expressions(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> None:
    """Plot the expression of the marker genes.

    Args:
        adata: annotated data matrix
        params: parameters

    """
    markers = pd.read_csv(params["markers"], sep="\t", index_col=False)
    marker_genes = list(set(markers["Markers"]).intersection(adata.var_names))

    marker_dict: dict[str, list[str]] = {}
    for k, v in dict(zip(markers["Markers"], markers["CellType"])).items():
        if k in marker_genes:
            marker_dict[v] = marker_dict.get(v, []) + [k]

    sc.settings.set_figure_params(dpi=300, figsize=[12, 8], fontsize=15)
    sc.pl.umap(adata, color=marker_genes, show=False)
    plt.savefig(params["marker_filename1"])
    plt.close()

    fig, axes = plt.subplots(3, 1, figsize=(15, 18))
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.2, hspace=0.3)
    for idx, key in enumerate(params["keys"]):
        sc.pl.dotplot(adata, marker_dict, key, dendrogram=True, show=False, ax=axes[idx])
    plt.savefig(params["marker_filename2"])
    plt.close()

    fig, axes = plt.subplots(3, 1, figsize=(15, 18))
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.2, hspace=0.3)
    for idx, key in enumerate(params["keys"]):
        sc.pl.matrixplot(
            adata,
            marker_dict,
            key,
            dendrogram=True,
            show=False,
            cmap=params["color_map"],
            standard_scale="var",
            colorbar_title="column scaled\nexpression",
            ax=axes[idx],
        )
    plt.savefig(params["marker_filename3"])
    plt.close()

    sc.pl.tracksplot(
        adata,
        marker_dict,
        groupby="leiden_1.0",
        dendrogram=False,
        show=False,
        figsize=(15, 8),
    )
    plt.savefig(params["marker_filename4"])
    plt.close()


def write_scdata(
    adata: ad.AnnData,
    h5ad: PathLike,
    h5seurat: PathLike,
) -> None:
    """Write scanpy analyzed data to H5AD and H5Seurat.

    Args:
        data: aggregated counts data
        h5ad: name of the h5ad
        h5seurat: name of the h5seurat data

    """
    pass


def plot_qc(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> None:
    """Plot QCs.

    Args:
        adata: annotated data matrix
        params: parameters

    """
    sc.pl.violin(
        adata,
        params["violin_plots"]["keys"],
        use_raw=params["violin_plots"]["use_raw"],
        groupby=params["violin_plots"]["groupby"],
        jitter=0.4,
        rotation=30,
        size=2,
        multi_panel=True,
        show=False,
    )
    plt.savefig(params["violin_plots"]["filename"])
    plt.close()

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    plt.subplots_adjust(left=0.1, right=0.7, bottom=0.1, top=0.9, wspace=0.2, hspace=0.3)
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        color="sample",
        show=False,
        use_raw=params["scatter_plots"]["use_raw"],
        ax=axes[0],
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="sample",
        show=False,
        use_raw=params["scatter_plots"]["use_raw"],
        ax=axes[1],
    )
    plt.savefig(params["scatter_plots"]["filename"])
    plt.close()

    sc.pl.highest_expr_genes(
        adata,
        show=False,
        n_top=params["box_plots"]["num_top_genes"],
    )
    plt.savefig(params["box_plots"]["filename"])
    plt.close()


def plot_sexbias(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> None:
    """Plot sex bias.

    Determine the sex of the sample by looking at reads from chrY (males)
    and XIST (X-inactive specific transcript) expression.

    Args:
        adata: annotated data matrix
        params: parameters

    """
    if set(params["sexbias_plots"]["keys"]) - set(adata.obs.columns):
        return

    sc.settings.set_figure_params(dpi=200, figsize=[8, 8], fontsize=8)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.2, hspace=0.3)
    for idx, key in enumerate(params["sexbias_plots"]["keys"]):
        sc.pl.violin(
            adata,
            key,
            use_raw=params["sexbias_plots"]["use_raw"],
            groupby=params["sexbias_plots"]["groupby"],
            jitter=0.4,
            rotation=30,
            size=2,
            multi_panel=True,
            show=False,
            ax=axes[idx],
        )
    plt.savefig(params["sexbias_plots"]["filename1"])
    plt.close()

    sc.settings.set_figure_params(dpi=200, figsize=[4, 4], fontsize=8)
    sc.pl.scatter(
        adata,
        x="XIST-counts",
        y="percent_chrY",
        color="sample",
        use_raw=params["sexbias_plots"]["use_raw"],
        show=False,
    )
    plt.savefig(params["sexbias_plots"]["filename2"])
    plt.close()


def correct_batch(
    adata: ad.AnnData,
    key: str = "runid",
) -> ad.AnnData:
    """Correct batches.

    Args:
        adata: annotated data matrix
        key: batch identifier

    Returns:
        batch-corrected annotated data

    """
    # create a new object with lognormalized counts
    adata_combat = sc.AnnData(X=adata.raw.X, var=adata.raw.var, obs=adata.obs)

    # first store the raw data
    adata_combat.raw = adata_combat
    sc.pp.combat(adata_combat, key=key)

    return adata_combat


def check_doublets(
    adata: ad.AnnData,
    params: ImmutableStrMapping,
) -> ad.AnnData:
    """Check doublets.

    Args:
        adata: annotated data matrix
        params: parameters

    Returns:
        annotated data with doublets information

    """
    scrub = scr.Scrublet(adata.X)
    adata.obs["doublet_scores"], adata.obs["predicted_doublets"] = scrub.scrub_doublets()
    scrub.plot_histogram()

    plt.savefig(params["filename"])
    plt.close()

    num_predicted = sum(adata.obs["predicted_doublets"])
    logging.info(f" ====> The number of predicted_doublets: {num_predicted}")
    logging.info(f' ====> The predicted_doublets_rate: {num_predicted / len(adata.obs["predicted_doublets"])}')

    return adata


def get_cluster_proportions(
    adata: ad.AnnData,
    cluster_key: str = "leiden_1.0",
    sample_key: str = "sample",
    exclude_samples: list[str] | None = None,
) -> pd.DataFrame:
    """Get cluster proportions.

    Args:
        adata: annotated data matrix
        cluster_key: adata.obs name storing clustering information
        sample_key: adata.obs name from adata, storing sample information
        exclude_samples: list of sample_key you like to exclude

    Returns:
        dataFrame with samples as the index and cluster proportion as values

    """
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=sample_key, index=cluster_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)

    if exclude_samples is not None:
        for sample in exclude_samples:
            props.drop(sample, axis=0, inplace=True)

    return props


def plot_cluster_proportions(
    cluster_props: pd.DataFrame,
    filename: PathLike,
    cluster_palette: str | None = None,
    xlabel_rotation: int = 90,
) -> None:
    """Plot cluster proportions.

    From Stacked barplot of scRNA-seq cluster proportions per sample:
    https://gist.github.com/wflynny/79c5266cc39a4a884958d696f84f85df

    """
    fig, ax = plt.subplots(dpi=300)
    fig.patch.set_facecolor("white")

    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette,
            n_colors=len(cluster_palette),
            as_cmap=True,
        )

    cluster_props.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        legend=None,
        colormap=cmap,
    )

    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title="Cluster Proportion")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    fig.tight_layout()

    plt.savefig(filename)
    plt.close()
