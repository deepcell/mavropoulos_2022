"""summarize.py for aggregating and summarizing the metrics and results."""

from __future__ import annotations

import argparse

from collections import defaultdict
from typing import Any, Callable, Mapping

import numpy as np
import pandas as pd

from scipy.stats import median_abs_deviation
from sklearn.decomposition import PCA

from .utils.pathutils import PathLike, SomePathLikes, parse_path, parse_paths
from .utils.plotutils import plot_boxplot, plot_histogram, plot_points


# FIXME: Just dtypes for now
# This will expand to include more, such as renaming, converter, and etc.
METADATA_DTYPES: dict[str, Any] = {
    "SampleName": str,
    "SequencingRun": str,
    "Type": str,
    "Well": str,
    "PreSortSpikeinRatio": np.float64,
    "LiveOrFixed": "category",
    "InputCellNum": np.float64,
    "InputDNApg": np.float64,
    "Panel": str,
    "Instrument": str,
    "SpikedinSamples": str,
    "Background": str,
    "JIRA": str,
    "Note": str,
}

METRICS_DTYPES: dict[str, Any] = {
    "sample": str,
    "n_frags": np.uint64,
    "median_frag": np.float64,
    "mad_frag": np.float64,
    "median_proper_frag": np.float64,
    "mad_proper_frag": np.float64,
    "pct_proper": np.float64,
    "pct_nonproper_same_chr": np.float64,
    "pct_nonproper_diff_chr": np.float64,
    "n_on_target_frags": np.uint64,
    "pct_on_target_frags": np.float64,
    "pct_on_target_bases": np.float64,
    "pct_ave_overlap": np.float64,
    "dup_rate": np.float64,
    "median_target_cov": np.float64,
    "q90_q10_ratio": np.float64,
    "pct_lt_5": np.float64,
}

COUNTS_DTYPES: dict[str, Any] = {
    "gene_id": str,
    "unstrand_cnt": np.uint64,
    "fwd_cnt": np.uint64,
    "rev_cnt": np.uint64,
}

RNA_METRICS_DTYPES: dict[str, Any] = {
    "SampleName": str,
    "NumFrags": np.uint64,
    "MismatchRate": np.float64,
    "MeanFragLen": np.float64,
    "MeanAdapterLen": np.float64,
    "PCTAdapter": np.float64,
    "PCT_A": np.float64,
    "PCT_G": np.float64,
    "PCT_UniqMapped": np.float64,
    "PCT_MultiMapped": np.float64,
    "PCT_Unmapped": np.float64,
    "PCT_FwdStranded": np.float64,
    "PCT_RevStranded": np.float64,
    "PCT_rRNA": np.float64,
}

FRAGMENT_DTYPES: dict[str, Any] = {
    "seqnames": str,
    "start": np.uint64,
    "end": np.uint64,
    "frag_len": np.uint64,
    "aln_type": "category",
    "amplicon": str,
    "on_or_off": "category",
    "on_target_bases": np.uint64,
    "dup_or_not": "category",
    "sample": str,
}

COVERAGE_DTYPES: dict[str, Any] = {
    "seqnames": str,
    "start": np.uint64,
    "end": np.uint64,
    "amplicon": str,
    "length": np.uint64,
    "gc": np.float64,
    "coverage": np.uint64,
    "norm_coverage": np.float64,
}

VARIANTS_DTYPES: dict[str, Any] = {
    "CHR": str,
    "POSITION": np.uint64,
    "REF": str,
    "ALT": str,
    "DEPTH": np.uint64,
    "RD": np.uint64,
    "AD": np.uint64,
    "AF": np.float64,
    "Sample": str,
    "Site": str,
    "VariantID": str,
}

VARIANTS_STAT_DTYPES: dict[str, Any] = {
    "MedianSiteCov": np.float64,
    "MinSiteCov": np.uint64,
    "MaxSiteCov": np.uint64,
    "NumSites": np.uint64,
    "NumALTCalls": np.uint64,
    "NumDropOut": np.uint64,
    "NumLT50Cov": np.uint64,
}

GENCODE_DTYPES: dict[str, Any] = {
    "chr": str,
    "start": np.uint64,
    "end": np.uint64,
    "info": str,
    "length": np.uint64,
    "strand": str,
    "gene_id": str,
    "gene_type": str,
    "gene_name": str,
    "anno_type": str,
}


def parse_to_list_string(value: str) -> list[str]:
    """Parse comma-delimited string to list of string.

    Args:
        value: the string value

    Returns:
        list of strings

    """
    return value.split(",")


def parse_to_list_int(value: str) -> list[int]:
    """Parse comma-delimited string to list of int values.

    Args:
        value: the string value

    Returns:
        list of ints

    """
    return [int(s) for s in value.split(",")]


def parse_to_list_float(value: str) -> list[float]:
    """Parse comma-delimited string to list of float values.

    Args:
        value: the string value

    Returns:
        list of floats

    """
    return [float(s) for s in value.split(",")]


def parse_percent(value: Any) -> float:
    """Parse percent into fraction.

    If **%** sign is present, remove **%** sign at the end in numerical value and return fraction.
    If **%** sign is not present, assume that **%** sign is not part of the value and just return fraction.

    Args:
        value: value to be parsed tin to fraction

    Returns:
        value converted into fraction

    Examples:
        >>> parse_percent("1 %")
        0.01
        >>> parse_percent("1%")
        0.01
        >>> parse_percent(1)
        0.01
        >>> parse_percent(0)
        0.0
        >>> parse_percent(np.nan)
        nan
        >>> parse_percent(np.inf)
        inf
        >>> parse_percent("NA")
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float: 'NA'
        >>> parse_percent("")
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float: ''
        >>> parse_percent(None)
        Traceback (most recent call last):
        ...
        TypeError: float() argument must be a string or a number, not 'NoneType'

    """
    if isinstance(value, str) and value.endswith("%"):
        value = value.rstrip("%").strip()

    return float(value) / 100


def parse_alignment_type_to_string(value: Any) -> str:
    """Parse alignment type numerical key to interpretable string output.

    The available mappings include:
        1 for Proper
        0 for non-proper on same chromosome
        -1 for non-proper on different chromosome

    Args:
        value representing 1, 0, or -1

    Returns:
        interpretable alignment type converted from numerical key

    Examples:
        >>> parse_alignment_type_to_string("1")
        'Proper'
        >>> parse_alignment_type_to_string(1.0)
        'Proper'
        >>> parse_alignment_type_to_string(1)
        'Proper'
        >>> parse_alignment_type_to_string("0")
        'NonProperSameChr'
        >>> parse_alignment_type_to_string(0.0)
        'NonProperSameChr'
        >>> parse_alignment_type_to_string(0)
        'NonProperSameChr'
        >>> parse_alignment_type_to_string("-1")
        'NonProperDiffChr'
        >>> parse_alignment_type_to_string(-1.0)
        'NonProperDiffChr'
        >>> parse_alignment_type_to_string(-1)
        'NonProperDiffChr'
        >>> parse_alignment_type_to_string(1.1)
        Traceback (most recent call last):
        ...
        KeyError: 1.1
        >>> parse_alignment_type_to_string(np.nan)
        Traceback (most recent call last):
        ...
        ValueError: cannot convert float NaN to integer
        >>> parse_alignment_type_to_string(None)
        Traceback (most recent call last):
        ...
        TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType'

    """
    if not isinstance(value, str) and value == int(value):
        value = str(int(value))
    return {
        "1": "Proper",
        "0": "NonProperSameChr",
        "-1": "NonProperDiffChr",
    }[value]


def get_highest_allele_frequency_variants(data: pd.DataFrame) -> pd.DataFrame:
    """Get highest allele frequency variants, if multiple variants exist.

    Args:
        variants data with allele frequency

    Returns:
        One variant per site, keep only highest AF, retain same behavior as seen from summarize.R

    """
    data = data.loc[~data["ALT"].isna()]
    return data[data.groupby([data.index, "Site"])["AF"].transform(max) == data["AF"]]


def compute_variants_stat(
    data: pd.DataFrame,
    no_depth_threshold: int = 0,
    low_depth_threshold: int = 50,
) -> pd.DataFrame:
    """Compute the variants summery stat.

    Args:
        data: aggregated variants data
        no_depth_threshold: depth at which is considered as drop out
        low_depth_threshold: depth at which is considered as low depth

    Returns:
        site coverage stats, number of call sites, no-ref calls, dropout sites, low-depth sites

    """
    zero_depth_filter = data["DEPTH"] == no_depth_threshold
    low_depth_filter = data["DEPTH"] <= low_depth_threshold
    no_alt_filter = data["ALT"] == "."

    coverage_stats = data.groupby(data.index).agg({"DEPTH": ["median", "min", "max"]})
    coverage_stats.columns = ["MedianSiteCov", "MinSiteCov", "MaxSiteCov"]

    n_sites = data.groupby(data.index)["Site"].nunique()
    n_nonref_calls = data.loc[~no_alt_filter].groupby(level=0)["Site"].count()
    n_snp_cov_dropout = data.loc[zero_depth_filter].groupby(level=0)["Site"].count()
    n_snp_lt_low_cov = data.loc[~zero_depth_filter & low_depth_filter].groupby(level=0)["Site"].count()

    n_stats = pd.concat([n_sites, n_nonref_calls, n_snp_cov_dropout, n_snp_lt_low_cov], axis=1)
    n_stats.columns = ["NumSites", "NumALTCalls", "NumDropOut", f"NumLT{low_depth_threshold}Cov"]

    return pd.concat([coverage_stats, n_stats], axis=1)


def compute_contamination_stat(
    data: pd.DataFrame,
    low_depth_threshold: int = 50,
) -> pd.DataFrame:
    """Compute contamination stat.

    From homozygous alt and homozygous ref sites, the observed minor allele count = contaminated
    FIXME: Not exactly the definition of contamination that we should be using. Need better one.

    Args:
        data: genotype with alt allele frequency data

    Returns:
        calculated contamination statistics

    """
    # FIXME: Remove fixme, when genotype inference is made. Currently not called.
    hom = data.loc[data["InferredGenotype"].isin(["0", "1"]) & (data["DEPTH"] > low_depth_threshold)]
    is_hom_alt = hom["InferredGenotype"] == "1"
    hom.loc[is_hom_alt, "AF"] = 1 - hom.loc[is_hom_alt, "AF"]

    return pd.DataFrame(
        {
            "ContamFrac": 2 * hom.groupby(hom.index)["AF"].mean(),
            "NumContamFracGT1pct": hom.loc[hom["AF"] > 0.01].groupby(level=0).size(),
            "NumContamFracGT10pct": hom.loc[hom["AF"] > 0.1].groupby(level=0).size(),
        }
    ).fillna(0)


def run_pca(
    data: pd.DataFrame,
    n_components: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run PCA analysis.

    Args:
        data: data to perform pca
        n_components: number of components used

    Returns:
        PC and variance

    """
    pca = PCA(n_components=n_components)
    pca_names = [f"PC{s+1}" for s in range(n_components)]
    components = pd.DataFrame(
        pca.fit_transform(data),
        columns=pca_names,
        index=data.index,
    )
    variance = pd.DataFrame(
        {"Variance": pca.explained_variance_ratio_},
        index=pca_names,
    )

    return components, variance


def compute_count_stat(counts: pd.DataFrame, metrics: pd.DataFrame) -> None:
    """Compute count stat.

    Args:
        counts: raw counts data to be updated with normalized counts
        metrics: metrics to be updated with mitochondrial representation and number of genes

    """
    norm_cnt_by_len = counts["fwd_cnt"] / counts["length"]
    size_factor = norm_cnt_by_len.groupby("sample").sum() / 1000000
    counts["norm_fwd_cnt"] = norm_cnt_by_len / size_factor

    chrM = counts.loc[counts["chr"] == "chrM"].reset_index().set_index("sample")
    metrics["PCT_mtRNA"] = chrM.groupby(chrM.index)["fwd_cnt"].agg("sum") / metrics["NumFrags"] * 100
    metrics["NumGeneDetected"] = counts.loc[counts["norm_fwd_cnt"] >= 0.1].groupby("sample").size()


def make_plots(
    output_dir: PathLike,
    *,
    agg_frags: pd.DataFrame,
    agg_covs: pd.DataFrame,
    pca_result: pd.DataFrame,
    variants: pd.DataFrame,
) -> None:
    """Make plots for fragments and coverage.

    Args:
        output_dir: directory for plots
        agg_frags: aggregated fragments data
        agg_covs: aggregated coverage data
        pca_result: PC1 and PC2 and % of variabilities explained by these
        variants: aggregated variants data

    """
    agg_frags = agg_frags.reset_index()
    agg_covs = agg_covs.reset_index()
    pca_result = pca_result.reset_index()

    output_dir = parse_path(output_dir)

    plot_histogram(
        agg_frags.loc[agg_frags["frag_len"] <= 5000],
        x="frag_len",
        fill="aln_type",
        col="sample",
        ncol=3,
        xlabel="Fragment Length",
    ).save(output_dir / "aggregated.frags.hist.pdf", height=4, width=8)

    plot_boxplot(
        agg_frags.loc[agg_frags["frag_len"] <= 5000],
        x="sample",
        y="frag_len",
        fill="aln_type",
        xlabel_rotation=65,
        ylabel="Fragment Length",
    ).save(output_dir / "aggregated.frags.box.pdf", height=4, width=8)

    plot_points(
        agg_covs,
        x="gc",
        y="coverage",
        col="sample",
        color="sample",
        ncol=3,
    ).save(output_dir / "aggregated.cov.gc_by_coverage.pdf", height=4, width=8)

    plot_points(
        agg_covs,
        x="gc",
        y="norm_coverage",
        col="sample",
        color="sample",
        ncol=3,
    ).save(output_dir / "aggregated.cov.gc_by_normcoverage.pdf", height=4, width=8)

    plot_points(
        agg_covs,
        x="length",
        y="coverage",
        col="sample",
        color="sample",
        ncol=3,
    ).save(output_dir / "aggregated.cov.len_by_coverage.pdf", height=4, width=8)

    plot_points(
        agg_covs,
        x="amplicon",
        y="norm_coverage",
        col="sample",
        color="sample",
        blank_x=True,
    ).save(output_dir / "aggregated.cov.targets.pdf", height=4, width=8)

    pca_result = pca_result.rename(columns={"index": "SampleName"})
    is_sample = pca_result["SampleName"] != "Variance"
    plot_points(
        pca_result.loc[is_sample],
        x="PC1",
        y="PC2",
        color="SampleName",
        xlabel=f"PC1: {round(pca_result.loc[~is_sample, 'PC1']*100, 2)}%",
        ylabel=f"PC2: {round(pca_result.loc[~is_sample, 'PC2']*100, 2)}%",
    ).save(output_dir / "aggregated.cov.pca.pdf")

    # FIXME: Uncomment when genotype calling is done
    # variants["AF2"] = np.where(variants["InferredGenotype"] == "1", 1 - variants["AF"], variants["AF"])
    # is_homozygous = variants["InferredGenotype"].isin(["0", "1"])
    # is_high_error = variants["AF2"] >= 0.01
    # variants.loc[is_homozygous & is_high_error, "InferredGenotype"] = "Contam"
    # plot_points(
    #     variants.loc[variants["InferredGenotype"] != "NA"].reset_index(),
    #     x="Ploidy",
    #     y="AF",
    #     color="InferredGenotype",
    #     xlabel="SNPs Ploidy",
    #     ylabel="Copy Number State",
    #     col="Sample",
    #     ncol=3,
    #     position="jitter",
    # ).save(output_dir / "genotype.pdf", width=8, height=4)


def aggregate_metrics(
    filenames: SomePathLikes,
    mapping: Mapping[str, Any],
    index_col: str = "sample",
    converters: dict[str, Callable[[str], Any]] | None = None,
    sampling_size: int | None = None,
) -> pd.DataFrame:
    """Aggregate sample metrics."""
    data = []
    for filename in sorted(parse_paths(filenames)):
        data.append(
            pd.read_csv(
                filename,
                dtype=mapping,
                index_col=index_col,
                converters=converters,
                low_memory=False,
            )
        )
    agg_data = pd.concat(data)

    # FIXME: Currently uses 100k frags, not all of them, sampling results in non-deterministic results
    if sampling_size is not None:
        agg_data = agg_data.sample(n=sampling_size)

    return agg_data


def aggregate_targets(
    filenames: SomePathLikes,
    mapping: Mapping[str, Any],
) -> pd.DataFrame:
    """Aggregate targets data."""
    data = []
    for filename in sorted(parse_paths(filenames)):
        samplename = filename.parent.name
        target = pd.read_csv(
            filename,
            dtype=mapping,
            sep="\t",
            index_col=None,
            low_memory=False,
        )
        target["sample"] = samplename
        data.append(target.reset_index(drop=True).set_index("sample"))

    return pd.concat(data)


def aggregate_counts(
    filenames: SomePathLikes,
    mapping: Mapping[str, Any],
) -> pd.DataFrame:
    """Aggregate counts data."""
    data = []
    for filename in sorted(parse_paths(filenames)):
        samplename = filename.parent.name
        counts = pd.read_csv(
            filename,
            dtype=mapping,
            index_col=0,
            sep="\t",
            header=None,
            names=mapping,
            low_memory=False,
        )
        counts["sample"] = samplename
        data.append(counts)

    return pd.concat(data)


def aggregate_variants(
    filenames: SomePathLikes,
    mapping: Mapping[str, Any],
) -> pd.DataFrame:
    """Aggregate variants data."""
    data = []
    for filename in sorted(parse_paths(filenames)):
        data.append(
            pd.read_csv(
                filename,
                dtype=mapping,
                sep="\t",
                index_col="Sample",
                converters={
                    "ALT": parse_to_list_string,
                    "AD": parse_to_list_int,
                    "AF": parse_to_list_float,
                },
                low_memory=False,
            )
        )
    agg_data = pd.concat(data)

    agg_data = agg_data.explode(["ALT", "AD", "AF"])
    agg_data["Site"] = agg_data["CHR"] + ":" + agg_data["POSITION"].astype(str) + ":" + agg_data["REF"]
    agg_data["VariantID"] = agg_data["Site"] + ":" + agg_data["ALT"]

    agg_data_high_freq = get_highest_allele_frequency_variants(agg_data)
    agg_variants_stat = compute_variants_stat(agg_data_high_freq)

    # FIXME: Do the genotype calling, mixture analysis during sample analysis in separate PR.
    # TODO: Add genotype to the calls
    # agg_data = add_genotype_samples(agg_data)
    # contamination = compute_contamination_stat(agg_data)
    # agg_variants_stat = pd.concat([agg_variants_stat, contamination], axis=1)

    return agg_data_high_freq, agg_variants_stat


def summarize_dna(args: argparse.Namespace, metadata: pd.DataFrame) -> None:
    """Summarize metrics."""
    # TODO: add metadata check, if sample name already exists or not.
    input_dir = args.input_dir / "metrics"

    metrics_dir = args.output_dir / "metrics"
    plots_dir = args.output_dir / "plots"

    metrics_dir.mkdir(exist_ok=True, parents=True)
    plots_dir.mkdir(exist_ok=True, parents=True)

    # read metrics and counts & aggregate
    agg_metrics = aggregate_metrics(input_dir.glob("**/*.metrics.csv"), METRICS_DTYPES)
    agg_fragments = aggregate_metrics(
        input_dir.glob("**/*frags.csv"),
        FRAGMENT_DTYPES,
        converters={"aln_type": parse_alignment_type_to_string},
        sampling_size=100000,
    )
    agg_coverage = aggregate_targets(input_dir.glob("**/*targets.tsv"), COVERAGE_DTYPES)
    agg_coverage_wide = agg_coverage.reset_index().pivot(
        index="amplicon",
        columns="sample",
        values="norm_coverage",
    )

    var_filenames = (args.input_dir / "germlines").glob("**/*snp.tsv")
    if args.workflow == "tgs":
        var_filenames = (args.input_dir / "mutations").glob("**/*mutation.tsv")
    agg_variants, agg_variants_stats = aggregate_variants(var_filenames, VARIANTS_DTYPES)

    pca_result, pca_variance = run_pca(agg_coverage_wide.T.fillna(1e-6))

    # Plot
    make_plots(
        plots_dir,
        agg_frags=agg_fragments,
        agg_covs=agg_coverage,
        pca_result=pd.concat([pca_result, pca_variance.T]),
        variants=agg_variants,
    )

    # Write target coverage and metrics summary
    agg_coverage.to_csv(metrics_dir / "aggregate.targets.cov.long.csv")
    agg_coverage_wide.to_csv(metrics_dir / "aggregate.targets.cov.wide.csv")
    agg_variants.to_csv(metrics_dir / "variants.long.tsv", sep="\t")

    pd.concat(
        [
            metadata,
            agg_metrics,
            agg_variants_stats,
        ],
        axis=1,
    ).to_csv(metrics_dir / "summary.tsv", index_label="SampleName", sep="\t")


def summarize_rna(args: argparse.Namespace, metadata: pd.DataFrame) -> None:
    """Summarize rna workflow metrics."""
    # FIXME: move to config
    input_dir = args.input_dir / "metrics"
    alignments_dir = args.input_dir / "alignments"

    metrics_dir = args.output_dir / "metrics"
    expression_dir = args.output_dir / "expression"

    metrics_dir.mkdir(exist_ok=True, parents=True)
    expression_dir.mkdir(exist_ok=True, parents=True)

    # read gene-level information
    expand_info = ["gene_id", "gene_type", "gene_name", "anno_type"]
    genes = pd.read_csv(
        args.genes,
        dtype=GENCODE_DTYPES,
        sep="\t",
        names=[s for s in GENCODE_DTYPES if s in set(GENCODE_DTYPES.keys()) - set(expand_info)],
        header=None,
        index_col=None,
        converters={"info": parse_to_list_string},
        low_memory=False,
    )

    genes[expand_info] = genes["info"].tolist()

    # read metrics and counts & aggregate
    agg_metrics = aggregate_metrics(
        input_dir.glob("**/*.metrics.csv"),
        RNA_METRICS_DTYPES,
        index_col="SampleName",
        converters={"MismatchRate": parse_percent},
    )

    agg_counts = aggregate_counts(alignments_dir.glob("**/*ReadsPerGene.out.tab"), COUNTS_DTYPES)
    # ignore N counts information, only use counts
    agg_counts = agg_counts.loc[~agg_counts.index.str.startswith("N_")]
    agg_counts = pd.merge(agg_counts, genes, on="gene_id").set_index(["sample", "gene_id"])

    compute_count_stat(agg_counts, agg_metrics)

    # write counts and metrics summary
    agg_counts.reset_index().pivot(
        index="gene_id",
        values="fwd_cnt",
        columns="sample",
    ).to_csv(expression_dir / "aggregate.gene.counts.tsv", sep="\t")

    pd.concat([metadata, agg_metrics], axis=1).to_csv(metrics_dir / "summary.tsv", sep="\t")


def get_wgs_cov_bias(
    base_line_coverage: pd.DataFrame,
    normalization_samples: list[str] | None = None,
) -> pd.DataFrame:
    """Calcualtes (chrom, bin) coverage biases using base line samples.

    Args:
        base_line_coverage: WGS14 normalization file.
        normalization_samples: list of normalization sample names.
        norm_header: header of the normalization file.

    Returns:
        Base stats for every (chrom, bin), including median, mad, and mad/median,
        calculated from all normalization samples.

    """
    norm_header = ["cov", "chrom", "pos_mb", "sample_name"]

    if normalization_samples is None:
        normalization_samples = [
            "GM12878-100c-Res-KAPA_S3",
            "GM12878-200c-Res-KAPA_S4",
            "GM24143-100c-Res-KAPA_S7",
            "GM24143-200c-Res-KAPA_S8",
            "GM24143-50c-Res-KAPA_S6",
            "GM24149-100c-Res-KAPA_S11",
            "GM24149-200c-Res-KAPA_S12",
            "GM24149-50c-Res-KAPA_S10",
            "GM24385-100c-Res-KAPA_S15",
            "GM24385-200c-Res-KAPA_S16",
            "GM24385-50c-Res-KAPA_S14",
            "HFL1-100c-Res-KAPA_S23",
            "HFL1-200c-Res-KAPA_S24",
            "HFL1-50c-Res-KAPA_S22",
        ]

    base_line_coverage.set_axis(norm_header, axis=1, inplace=True)

    base_line_coverage = get_wgs_scaled_coverage(
        base_line_coverage,
        samples=normalization_samples,
    )

    base_line_coverage.dropna(inplace=True)
    base_line_coverage.reset_index(inplace=True)
    base_line_coverage.drop(columns=["sample_mean"])

    uniques = base_line_coverage.iloc[
        base_line_coverage.duplicated(subset=["chrom", "pos_mb"])[:].index[
            ~base_line_coverage.duplicated(subset=["chrom", "pos_mb"])[:]
        ]
    ]

    result: defaultdict[str, list[Any]] = defaultdict(list)

    for _, row in uniques.iterrows():

        tmp = base_line_coverage.loc[
            (base_line_coverage["chrom"] == row["chrom"]) & (base_line_coverage["pos_mb"] == row["pos_mb"])
        ]

        mad = median_abs_deviation(tmp["scaled_cov"], scale="normal")

        result["chrom"].append(row["chrom"])
        result["pos_mb"].append(row["pos_mb"])
        result["bias_median"].append(tmp["scaled_cov"].median())
        result["bias_mad"].append(mad)
        result["bias_cv"].append(mad / tmp["scaled_cov"].median())

    return pd.DataFrame.from_dict(result)


def get_wgs_scaled_coverage(
    genome_mb_coverage: pd.DataFrame,
    samples: list[str] | None = None,
) -> pd.DataFrame:
    """Calculate chromosomal scaled coverage information.

    Args:
        genome_mb_coverage: Input coverage to be scaled.
        samples: list of samples to be normalized.

    Returns:
        Normalized samples.

    """
    header = ["cov", "chrom", "pos_mb", "sample_name"]

    genome_mb_coverage.set_axis(header, axis=1, inplace=True)
    sample_name = header[3]
    cov = header[0]

    if samples is None:
        samples = genome_mb_coverage[sample_name].unique()

    for sample in samples:
        sample_rows = genome_mb_coverage.loc[genome_mb_coverage[sample_name] == sample]
        sample_coverages = sample_rows[cov]
        coverage_mean = np.mean(sample_coverages)
        genome_mb_coverage.loc[genome_mb_coverage[sample_name] == sample, "sample_mean"] = coverage_mean
        genome_mb_coverage.loc[genome_mb_coverage[sample_name] == sample, "scaled_cov"] = (
            sample_coverages / coverage_mean
        )

    return genome_mb_coverage


def merge_bias_cov(
    cov_bias: pd.DataFrame,
    scaled_cov: pd.DataFrame,
    bias_cv_threshold: float = 0.2,
    include_chromosome: str | None = None,
) -> pd.DataFrame:
    """Merge coverage bias with scaled coverage by (chrom, pos_mb).

    Args:
        cov_bias: coverage bias from base normalization samples.
        scaled_cov: sample scaled coverages.
        bias_cv_threshold: bias cv threshold value, those above which will be ignored.
        include_chromosome: pattern of chromosomes to be included, include all chromosomes if pattern not given.

    Returns:
        Merged and numbered.

    """
    if include_chromosome is not None:
        cov_bias = cov_bias[cov_bias["chrom"].str.contains(include_chromosome, regex=True)]
        scaled_cov = scaled_cov[scaled_cov["chrom"].str.contains(include_chromosome, regex=True)]

    result = scaled_cov.merge(cov_bias, "inner", on=["chrom", "pos_mb"])

    if result.shape[0] == 0:
        raise ValueError("No data left in the merged dataframe after filtering.")

    result["chrom_num"] = result.apply(
        lambda row: int(row["chrom"].split("chr")[1]),
        axis=1,
    )
    result.drop(result[result["bias_cv"] >= bias_cv_threshold].index, inplace=True)

    return result


def summarize(args: argparse.Namespace) -> None:
    """Summarize metrics and variants."""
    if args.workflow == "rna" and args.genes is None:
        raise ValueError("Requires gencode annotation for RNA metrics summary.")

    metadata = pd.read_csv(args.metasheet, dtype=METADATA_DTYPES, sep="\t", index_col=0, low_memory=False)
    summarize_rna(args, metadata) if args.workflow == "rna" else summarize_dna(args, metadata)


def arg_parser() -> argparse.ArgumentParser:
    """Create an argument parser for summarize.py."""
    parser = argparse.ArgumentParser(description="Summarize metrics and results.")

    parser.add_argument("metasheet", metavar="FILE", type=parse_path, help="metasheet")
    parser.add_argument("input_dir", metavar="DIR", type=parse_path, help="input directory")
    parser.add_argument("-o", "--output_dir", metavar="DIR", required=True, type=parse_path, help="output directory")
    parser.add_argument("--workflow", metavar="STR", help="analysis workflow type: sid, tgs, or rna")
    parser.add_argument("--genes", metavar="TSV", type=parse_path, help="gencode annotation at gene level")
    parser.set_defaults(func=summarize)

    return parser


def main(
    *,
    argv: list[str] | None = None,
    args: argparse.Namespace | None = None,
) -> None:
    """Summarize CLI."""
    parser = arg_parser()

    if argv is not None and args is not None:
        raise ValueError("argv and args are mutually exclusive")
    elif args is None:
        args = parser.parse_args(argv)

    if args.func:
        args.func(args)
    else:
        parser.print_help()
