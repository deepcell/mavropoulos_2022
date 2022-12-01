"""wrappers.py wrapper functions to build command, API abstraction for CLI."""

from __future__ import annotations

from typing import Iterable

from .utils.pathutils import PathLike, SomePathLikes, parse_path, parse_paths
from .utils.shellutils import Command


def make_command_atropos_trim(
    *,
    algorithm: str,
    adapter1: str | Iterable[str],
    read1: PathLike,
    trimmed_read1: PathLike,
    adapter2: str | Iterable[str] | None = None,
    read2: PathLike | None = None,
    trimmed_read2: PathLike | None = None,
    trim_info: PathLike | None = None,
    trim_report: PathLike | None = None,
    logfile: PathLike | None = None,
    loglevel: str = "INFO",
    no_cache_adapters: bool = True,
    nextseq_trim_qscore: int | None = None,
    minimum_length: int | None = 50,
    maximum_length: int | None = None,
    quality_cutoff: int | None = 20,
    max_n: int | None = 5,
    trim_n: bool = True,
    op_order: str = "GAWCQ",
    sampleid: str | None = None,
    threads: int = 2,
) -> Command:
    """Create command for atropos, a fork of cutadapt 1.10.

    https://github.com/marcelm/cutadapt/tree/2f3cc0717aa9ff1e0326ea6bcb36b712950d4999

    Args:
        algorithm: algorithm to be used during the trimming, one of adapter or insert
        adapter1: adapter sequence to be removed from read1
        read1: demultiplexed read, either from single-end or 1st from paired-end read
        trimmed_read1: adapter1 removed read
        adapter2: adapter sequence to be removed from read2
        read2: demultiplexed read, 2nd from paired-end read
        trimmed_read2: adapter2 removed read
        trim_info: information about each read and its adapter matches
        trim_report: written report rather than stdout/stder
        logfile: file with logging information
        loglevel: logging level, one of DEBUG, INFO (default), WARN, ERROR
        no_cache_adapters: do not fache adapters list as `.adapters` in the working directory
        nextseq_trim_qscore: if set, perform nextseq-specific quality trim with dark cycle
        minimum_legnth: minimum allowed trimmed reads
        maximum_length: maximum allowed trimmed reads
        quality_cutoff: low-quality bases
        max_n: maximum allowed N bases before discarded
        trim_n: if set True, trim N's on the ends of reads
        op_order: order of trim operation
        sampleid: sample ID, if defined, this id is used for metrics summary output
        threads: number of threads to use for read trimming, 0 uses max available

    Returns:
        trimming command

    """
    cmd = Command(
        "atropos",
        "trim",
        "-T",
        threads,
        "--aligner",
        algorithm,
        "--op-order",
        op_order,
        "--preserve-order",
    )

    if minimum_length is not None:
        cmd += ["--minimum-length", minimum_length]

    if maximum_length is not None:
        cmd += ["--maximum-length", maximum_length]

    if quality_cutoff is not None:
        cmd += ["--quality-cutoff", quality_cutoff]

    if max_n is not None:
        cmd += ["--max-n", max_n]

    if trim_n:
        cmd += ["--trim-n"]

    if no_cache_adapters:
        cmd += ["--no-cache-adapters"]

    if nextseq_trim_qscore is not None:
        cmd += ["--nextseq-trim", nextseq_trim_qscore]

    if isinstance(adapter1, str):
        adapter1 = [adapter1]

    if adapter2 is not None and isinstance(adapter2, str):
        adapter2 = [adapter2]

    if adapter2 is None:
        for a1 in adapter1:
            cmd += ["--adapter", a1]

        cmd += [
            "-se",
            read1,
            "--output",
            trimmed_read1,
        ]

    elif read2 is not None and trimmed_read2 is not None:
        for a1 in adapter1:
            cmd += ["--adapter", a1]

        for a2 in adapter2:
            cmd += ["--adapter2", a2]

        cmd += [
            "-pe1",
            read1,
            "-pe2",
            read2,
            "--output",
            trimmed_read1,
            "--paired-output",
            trimmed_read2,
        ]

    else:
        raise ValueError("Both fastq and trimmed output filenames are needed for paired end data")

    if sampleid is not None:
        cmd += ["--sample-id", sampleid]

    if trim_info is not None:
        cmd += ["--info-file", trim_info]

    if trim_report is not None:
        cmd += ["--report-file", trim_report]

    if logfile is not None:
        cmd += [
            "--log-file",
            logfile,
            "--log-level",
            loglevel,
        ]

    return cmd


def make_command_star_alignreads(
    *,
    reference: PathLike,
    read1: PathLike,
    aligned_read: PathLike,
    read2: PathLike | None = None,
    reference_spikein: PathLike | None = None,
    gzip_compressed: bool = True,
    quant_mode: str = "GeneCounts",
    threads: int = 1,
    tempdir: PathLike | None = None,
    outdir: PathLike | None = None,
) -> Command:
    """Create command for star alignment.

    STAR pre-determines the name of the output files, using the prefix given.
    The outputs will be renamed, once STAR completes  the alignment.

    Args:
        reference: directory where reference genome files are stored
        reference_spikein: reference genome for spike-in sample
        read1: read 1
        read2: read 2
        gzip_compressed: gzip compressed, TODO: support bzip2
        aligned_read: aligned read
        quant_mode: type of quantification, TranscriptomeSAM or GeneCounts (default)
        tempdir: if set, use this directory as temporary, else default to outFileNamePrefix_STARtmp
        outdir: output directory with log and alignment file

    Returns:
        STAR alignment command

    """
    aligned_read_default = f"{parse_path(aligned_read).parent}/"
    reads = [read1] if read2 is None else [read1, read2]

    cmd = Command(
        "STAR",
        "--runMode",
        "alignReads",
        "--runThreadN",
        threads,
        "--genomeDir",
        reference,
        "--quantMode",
        quant_mode,
        "--outFileNamePrefix",
        aligned_read_default,
        "--readFilesIn",
        *reads,
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        "--outReadsUnmapped",
        "FastX",
        "--outSAMunmapped",
        "Within",
        "--outSAMattributes",
        "All",
    )

    if reference_spikein is not None:
        cmd += ["--genomeFastaFiles", reference_spikein]

    if gzip_compressed:
        cmd += ["--readFilesCommand", "zcat"]

    return cmd


def make_command_bwa_mem(
    *,
    reference: PathLike,
    read1: PathLike,
    aligned_read: PathLike,
    read2: PathLike | None = None,
    readgroup: str | None = None,
    mark_secondary: bool = True,
    temp_dir: PathLike | None = None,
    min_mapq: int | None = None,
    process_n_bases: int = 1000000,
    threads: int = 1,
    logfile: PathLike | None = None,
) -> Command:
    r"""Make command for bwa mem alignment, sort and return.

    Process 1000000 input bases in each batch regardless of nThreads for reproducibility.

    FIXME: WARNING: With given command, the alignments are non-deterministic.

    Args:
        reference: reference genome
        read1: read1 fastq
        aligned_read: bwa mem aligned with parameters given
        read2: read2 fastq
        readgroup: read group header line, in the format of '@RG\tID:foo\tSM:bar'
        mark_secondary: if set True, mark shorter split hits as secondary
        min_mapq: only include reads with mapping quality >= INT
        process_n_bases: process N input bases in each batch for reproducibility, regardless of nThreads
        threads: number of threads
        logfile: file with logging information

    Returns:
        bwa alignment command

    """
    # tmp_read = parse_path(temp_dir) / "temp.bam"
    reads = [read1] if read2 is None else [read1, read2]

    cmd1 = Command("bwa", "mem", "-t", threads, "-K", process_n_bases, reference, *reads)

    if mark_secondary:
        cmd1 += ["-M"]

    if readgroup is not None:
        cmd1 += ["-R", readgroup]

    cmd1 = cmd1.redirect_output(stderr=logfile)
    cmd2 = make_command_samtools_view(min_mapq=min_mapq)
    cmd3 = make_command_samtools_sort(sorted_read=aligned_read, threads=threads)

    return cmd1 | cmd2 | cmd3


def make_command_samtools_sort(
    *,
    aligned_read: PathLike | None = None,
    sorted_read: PathLike | None = None,
    output_format: str = "BAM",
    temp_read: PathLike | None = None,
    sort_by_name: bool = False,
    threads: int = 1,
) -> Command:
    """Make command for samtools sort.

    Args:
        aligned_read: aligned read to be sorted
        sorted_read: sorted aligned read
        output_format: one of SAM, BAM (default), CRAM
        temp_read: temporary output to PREFIX.nnnn.bam
        sort_by_name: sort by read name, not by coordinate
        threads: number of additional theads to use

    Returns:
        samtools sort command

    """
    cmd = Command(
        "samtools",
        "sort",
        "--threads",
        threads,
    )

    if sorted_read is not None:
        cmd += ["-o", sorted_read, "--output-fmt", output_format]

    if temp_read is not None:
        cmd += ["-T", temp_read]

    if sort_by_name:
        cmd += ["-n"]

    if aligned_read is not None:
        cmd += [aligned_read]

    return cmd


def make_command_samtools_view(
    *,
    aligned_read: PathLike | None = None,
    output_read: PathLike | None = None,
    output_format: str = "BAM",
    reference: PathLike | None = None,
    include_header: bool = True,
    min_mapq: int | None = None,
    include_flags: int | None = None,
    exclude_flags: int | None = None,
) -> Command:
    """Make command for samtools view.

    Args:
        aligned_read: input aligned read in BAM,SAM, CRAM
        output_read: output file name
        output_format: one of BAM (default), SAM, or CRAM
        reference: reference FASTA required for CRAM output
        include_header: include header in the output
        min_mapq: only include reads with mapping quality >= INT
        include_flags: only include reads with all of the FLAGs in INT present
        exclude_flags: only include reads with none of the FLAGS in INT present

    Returns:
        samtools view command

    Raises:
        ValueError: if reference is not provided with CRAM output

    """
    if output_format == "CRAM" and reference is None:
        raise ValueError("CRAM requires reference sequence FASTA file")

    cmd = Command("samtools", "view", "-C" if output_format == "CRAM" else "-b")

    if include_header:
        cmd += ["-h"]

    # FIXME: is it wrong to specify reference for format!=CRAM?
    if output_format == "CRAM" and reference is not None:
        cmd += ["--reference", reference]

    if output_read is not None:
        cmd += [
            "-o",
            output_read,
            "--output-fmt",
            output_format,
        ]

    if min_mapq is not None:
        cmd += ["-q", min_mapq]

    if include_flags is not None:
        cmd += ["-f", include_flags]

    if exclude_flags is not None:
        cmd += ["-F", exclude_flags]

    if aligned_read is not None:
        cmd += [aligned_read]

    return cmd


def make_command_samtools_index(
    *,
    aligned_read: PathLike,
    format: str = "bai",
    threads: int = 1,
) -> Command:
    """Make command for samtools index.

    Args:
        aligned_read: aligned read filename
        format: BAI-format index for BAM, else CSI-format
        threads: number of of threads

    Returns:
        samtools indexing read command

    Raises:
        ValueError: requested format for index is neither bai nor csi

    """
    if format not in ("bai", "csi"):
        raise ValueError(f"{format=} is not supported. It can either be bai or csi.")

    return Command(
        "samtools",
        "index",
        "-@",
        threads,
        "-b" if format == "bai" else "-c",
        aligned_read,
    )


def make_command_samtools_flagstat(
    *,
    aligned_read: PathLike,
    stat: PathLike,
    format: str = "json",
    threads: int = 1,
) -> Command:
    """Make command for samtools flagstat.

    Args:
        aligned_read: aligned read filename
        stat: flag stat output
        format: json or tsv output format
        threads: number of of threads

    Returns:
        samtools flagstat command

    """
    if format not in ("tsv", "json"):
        raise ValueError(f"{format=} is not supported (tsv or json are allowed)")

    return Command(
        "samtools",
        "flagstat",
        "--threads",
        threads,
        "--output-fmt",
        format,
        aligned_read,
    ).redirect_output(stdout=stat)


def make_command_samtools_fixmate(
    *,
    aligned_read: PathLike,
    fixed_read: PathLike,
    threads: int = 1,
) -> Command:
    """Make command for samtools fixmate.

    Args:
        aligned_read: aligned read filename
        fixed_read: mate coordinates and insert size size fixed aligned read filename
        threads: number of of threads

    Returns:
        samtools fixmate command

    Raises:
        ValueError: requested format for index is neither bai nor csi

    """
    sort = make_command_samtools_sort(
        aligned_read=aligned_read,
        sort_by_name=True,
    )

    fixmate = Command(
        "samtools",
        "fixmate",
        "--threads",
        threads,
        "-",
        fixed_read,
    )

    return sort | fixmate


def make_command_samtools_calmd(
    *,
    aligned_read: PathLike,
    reference: PathLike,
    threads: int = 1,
) -> Command:
    """Make command for samtools calmd to fix NM (edit distance) & MD (mismatch and deleted bases) tags.

    Args:
        aligned_read: sorted, aligned reads
        reference: reference FASTA required for CRAM output
        threads: number of of threads

    Returns:
        samtools calmd command

    """
    return Command(
        "samtools",
        "calmd",
        "--threads",
        threads,
        "-b",
        aligned_read,
        reference,
    )


def make_command_samtools_ampliconclip(
    *,
    reference: PathLike,
    aligned_read: PathLike,
    clipped_read: PathLike,
    primer: PathLike,
    tmp_dir: PathLike,
    hard_clip: bool = False,
    both_ends: bool = True,
    fail_len: int = 50,
    threads: int = 1,
) -> Command:
    """Make command for samtools ampliconclip.

    Clipping reads results in template length (TLEN) being incorrect, as well as any MD and NM tags.
    samtools sort, fixmate, and calmd commands are run together to correct.

    Args:
        reference: reference genome
        aligned_read: aligned read filename
        clipped_read: amplicon clipped aligned read
        primer: primer target coordinate in BED
        tmp_dir: store intermediate bam files
        hard_clip: if set True, hard clip amplicon primers, default is soft clip
        both_ends: if set True, clip on both 5' and 3' ends, default is just on 5' ends
        fail_len: mark as QCFAIL if shorter than specified
        threads: number of threads

    Returns:
        samtools ampliconclip command

    """
    temp_dir = parse_path(tmp_dir)
    tmp_clipped = temp_dir / "tmp_clipped.bam"
    tmp_fixed = temp_dir / "tmp_fixed.bam"

    clip = Command(
        "samtools",
        "ampliconclip",
        "--threads",
        threads,
        "--fail-len",
        fail_len,
        "-b",
        primer,
        aligned_read,
        "-o",
        tmp_clipped,
    )

    if hard_clip:
        clip += ["--hard-clip"]

    if both_ends:
        clip += ["--both-ends"]

    fixmate = make_command_samtools_fixmate(
        aligned_read=tmp_clipped,
        fixed_read=tmp_fixed,
    )

    sort = make_command_samtools_sort(aligned_read=tmp_fixed)

    calmd = make_command_samtools_calmd(
        aligned_read="-",
        reference=reference,
    )

    return clip & fixmate & sort | calmd.redirect_output(stdout=clipped_read)


def make_command_filter_mutect2(
    reference: PathLike,
    variants: PathLike,
    filtered_calls: PathLike,
    stats: PathLike | None = None,
    min_unique_alt_cnt: int | None = None,
    min_reads_per_strand: int | None = None,
    max_events_in_region: int | None = None,
    min_allele_fraction: float | None = None,
) -> Command:
    """Make command for filtering mutect2 calls made.

    # FIXME: log this later

    Args:
        reference: reference sequence file
        variants: A VCF file containing variants
        filtered_calls: The output filtered VCF file
        stats: mutect stats file output by Mutect2
        min_unique_alt_cnt: min unique reads supporting the alternate alleles
        min_reads_per_strand: min alt reads required on both forward and reverse strands
        max_events_in_region: max events in a single assembly region
        min_allele_fraction: min allele fraction required

    Returns:
        gatk FilterMutectCalls command

    """
    cmd = Command(
        "gatk",
        "FilterMutectCalls",
        "--reference",
        reference,
        "--variant",
        variants,
        "--output",
        filtered_calls,
    )

    if min_unique_alt_cnt is not None:
        cmd += ["--unique-alt-read-count", min_unique_alt_cnt]

    if min_reads_per_strand is not None:
        cmd += ["--min-reads-per-strand", min_reads_per_strand]

    if max_events_in_region is not None:
        cmd += ["--max-events-in-region", max_events_in_region]

    if min_allele_fraction is not None:
        cmd += ["--min-allele-fraction", min_allele_fraction]

    if stats is not None:
        cmd += ["--stats", stats]

    return cmd


def make_command_call_mutect2(
    *,
    reference: PathLike,
    aligned_reads: SomePathLikes,
    calls: PathLike,
    intervals: PathLike | None = None,
    min_mapping_quality: int = 30,
    min_base_quality: int = 25,
    min_length: int = 50,
    max_length: int = 1000,
    max_reads_per_alignment_start: int = 0,
    filters: str | Iterable[str] | None = None,
    evidence_bam: PathLike | None = None,
    normal_sample_ids: str | Iterable[str] | None = None,
    tumor_sample_ids: str | Iterable[str] | None = None,
    germline_resource: PathLike | None = None,
    make_haplotype_graph: bool = False,
    dont_use_soft_clipped_bases: bool = False,
) -> Command:
    """Make command for GATK Mutect2.

    # FIXME: Log this

    Args:
        reference: path to reference sequence file
        aligned_reads: BAM/SAM/CRAM file containing reads. May be specified multiple times
        calls: output variant calls
        intervals: file with genomic intervals over which to operate
        min_mapping_quality: minimum mapping quality to keep (inclusive)
        min_base_quality: minimum base quality required to consider a base for calling
        min_length: keep only reads with length at least equal to the specified value
        max_length: keep only reads with length at most equal to the specified value
        max_reads_per_alignment_start: max# reads to retain per alignment start position. Use 0 to disable
        filters: read filters to be applied before analysis. May be specified multiple times
        evidence_bam: file to which assembled haplotypes should be written
        normal_sample_ids: BAM sample name of normal(s). May be specified multiple times
        tumor_sample_ids: BAM sample name of tumor sample(s). May be specified multiple times
        germline_resource: population vcf of germline sequencing containing allele fraction
        make_haplotype_graph: if set True, construct a Linked De Bruijn graph
        dont_use_soft_clipped_bases: if set True, don't analyze soft clipped bases in the reads

    Returns:
        Mutect2 command for paired normal/tumor variant call

    """
    cmd = Command(
        "gatk",
        "Mutect2",
        "--reference",
        reference,
        "--output",
        calls,
        "--minimum-mapping-quality",
        min_mapping_quality,
        "--min-base-quality-score",
        min_base_quality,
        "--min-read-length",
        min_length,
        "--max-read-length",
        max_length,
        "--max-reads-per-alignment-start",
        max_reads_per_alignment_start,
    )

    for aligned_read in parse_paths(aligned_reads):
        cmd += ["--input", aligned_read]

    if intervals is not None:
        cmd += ["--intervals", intervals]

    if isinstance(filters, str):
        filters = [filters]

    for read_filter in filters or []:
        cmd += ["--read-filter", read_filter]

    if dont_use_soft_clipped_bases:
        cmd += ["--dont-use-soft-clipped-bases"]

    if make_haplotype_graph:
        cmd += ["--linked-de-bruijn-graph"]

    if germline_resource is not None:
        cmd += ["--germline-resource", germline_resource]

    if evidence_bam is not None:
        cmd += ["--bamout", evidence_bam]

    if isinstance(normal_sample_ids, str):
        normal_sample_ids = [normal_sample_ids]

    if isinstance(tumor_sample_ids, str):
        tumor_sample_ids = [tumor_sample_ids]

    for normal_sample_id in normal_sample_ids or []:
        cmd += ["--normal-sample", normal_sample_id]

    for tumor_sample_id in tumor_sample_ids or []:
        cmd += ["--tumor-sample", tumor_sample_id]

    return cmd


def make_command_call_germline(
    *,
    reference: PathLike,
    aligned_read: PathLike,
    targets: PathLike,
    calls: PathLike,
    sid: str,
    save_pileup: bool = False,
    batchsize: int = 500,
    threads: int = 1,
) -> Command:
    """Make command for germline call in naive mode.

    # TODO: log this

    Args:
        reference: reference sequence
        aligned_read: SAM/BAM/CRAM file
        targets: the targeted regions (BED) or positions (VCF)
        calls: somatic mutations calls output filename
        sid: sample identifier
        save_pileup: if set true, output pileup file
        batchsize: the number of regions/positions to process in parallel
        threads: the number of threads

    Returns:
        naive caller command

    """
    cmd = Command(
        "call_naive",
        "--nproc",
        threads,
        "--bam",
        aligned_read,
        "--ref",
        reference,
        "--target",
        targets,
        "--sample",
        sid,
        "--out",
        calls,
        "--batch_size",
        batchsize,
    )

    if save_pileup:
        cmd += ["--pileup"]

    return cmd


def make_command_bcftools_mpileup(
    *,
    reference: PathLike,
    aligned_reads: SomePathLikes,
    target: PathLike | None = None,
    output: PathLike | None = None,
    output_type: str | None = None,
    annotate: str | None = None,
    include_flags: str | int | None = None,
    exclude_flags: str | int | None = None,
    max_depth: int = 999999,
    min_mapq: int = 30,
    min_baseq: int = 25,
    threads: int = 1,
) -> Command:
    """Make command for sample swap check.

    # TODO: log this

    Args:
        reference: faidx indexed reference sequence file
        aligned_reads: 1 or more aligned read input file(s)
        target: when used, restrict to regions listed in a file
        output: file output of pileup data
        output_type: BCF, 'u' or 'b' for compressed; VCF, 'v' (default) or 'z' for compressed
        include_flags: required flags, skip reads with mask bits unset
        exclude_flags: filter flags, skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]
        max_depth: max per-file depth; avoids excessive memory usage
        min_mapq: skip alignments if minimum mapQ requirement is not met
        min_baseq: skip bases if minimum baseQ/BAQ requirement is not met
        annotate: comma separated optional tags to output
        threads: number of extra output compression threads

    Returns:
        bcftools pileup command

    """
    cmd = Command(
        "bcftools",
        "mpileup",
        "--threads",
        threads,
        "--fasta-ref",
        reference,
        "--min-MQ",
        min_mapq,
        "--min-BQ",
        min_baseq,
        "--max-depth",
        max_depth,
    )

    if target is not None:
        cmd += [
            "--regions-file",
            target,
        ]

    if output is not None:
        cmd += ["--output", output]

    if annotate is not None:
        cmd += ["--annotate", annotate]

    if output_type is not None:
        if output_type not in ("u", "v", "b", "z"):
            raise ValueError(f"{output_type=} is not supported, must be one of b, u, v, z")
        cmd += ["--output-type", output_type]

    if include_flags is not None:
        cmd += ["--incl-flags", include_flags]

    if exclude_flags is not None:
        cmd += ["--excl-flags", exclude_flags]

    cmd += parse_paths(aligned_reads)

    return cmd


def make_command_bcftools_call(
    input: PathLike | None = None,
    output: PathLike | None = None,
    output_type: str | None = None,
    use_multiallelic: bool = False,
    out_variants_only: bool = False,
) -> Command:
    """Make command for bcftools call.

    # FIXME: log this

    Args:
        input: VCF or BCF input file to be used to make calls
        output: SNP/indel variant calls output
        output_type: BCF, 'u' or 'b' for compressed; VCF, 'v' (default) or 'z' for compressed
        use_multiallelic: alternative model for multiallelic and rare-variant calling
        out_variants_only: output variant sites only

    Returns:
        bcftools call command

    """
    cmd = Command("bcftools", "call")

    if output is not None:
        cmd += ["--output", output]

    if output_type is not None:
        cmd += ["--output-type", output_type]

    if use_multiallelic:
        cmd += ["--multiallelic-caller"]

    if out_variants_only:
        cmd += ["--variants-only"]

    if input is not None:
        cmd += [input]

    return cmd


def make_command_check_swap_metrics(
    input: PathLike,
    output: PathLike,
) -> Command:
    """Make command to check metrics for sample swap.

    Args:
        input: VCF input file from variant calls from bcftools pileup
        output: metrics for sample swapping

    Returns:
        command to generate sample swap metrics

    """
    return Command("check_swap", "--vcf", input, "--out", output)


def make_command_collect_metrics(
    *,
    input_type: str,
    reference: PathLike,
    aligned_read: PathLike,
    metrics_dir: PathLike,
    sid: str,
    target: PathLike | None = None,
    annotation: PathLike | None = None,
) -> Command:
    """Collect metrics.

    # FIXME: log this & change to just API calling to aggregate?

    Args:
        input_type: input material type, either DNA or RNA
        reference: reference genome in fasta
        aligned_read: aligned read file in BAM (TODO: expand to CRAM)
        metrics_dir: metrics output directory
        sid: sample identifier
        target: DNA panel target used, only required if input type is DNA
        annotation: gtf annotation file, only required if input type is RNA

    Returns:
        command to run metrics aggregation

    """
    if input_type == "DNA" and target is None:
        raise ValueError("DNA metrics aggregation requires the target information")

    if input_type == "RNA" and annotation is None:
        raise ValueError("RNA metrics aggregation requires the GTF annotation file")

    cmd = Command(
        "collect_metrics",
        input_type,
        "--bam",
        aligned_read,
        "--sample",
        sid,
        "--out",
        metrics_dir,
        "--genome",
        reference,
    )

    if target is not None:
        cmd += ["--panel", target]

    if annotation is not None:
        cmd += ["--gtf", annotation]

    return cmd


def make_command_bcftools_index(
    filename: PathLike,
    *,
    use_csi: bool = True,
    output: PathLike | None = None,
    threads: int = 1,
) -> Command:
    """Make command for indexing variant file.

    Args:
        filename: name of the VCF/BCF file to be indexed, <in.bcf>|<in.vcf.gz>
        use_csi: if set True (default), generate CSI index for VCF, intead of TBI
        output: optional output index file name
        threads: number of threads to be used

    Returns:
        bcftools index command

    """
    cmd = Command("bcftools", "index", "--threads", threads)

    if use_csi:
        cmd += ["--csi"]
    else:
        cmd += ["--tbi"]

    if output is not None:
        cmd += ["--output-file", output]

    cmd += [filename]

    return cmd


def make_command_tabix(
    *,
    filename: PathLike,
    use_csi: bool = False,
) -> Command:
    """Make command for indexing variant file.

    Args:
        filename: name of the file to be indexed
        use_csi: if set True, generate CSI index for VCF, intead of TBI

    Returns:
        tabix command

    """
    cmd = Command("tabix", filename)

    if use_csi:
        cmd += ["--csi"]

    return cmd


def build_command_shuffle(
    filename: PathLike,
    output: PathLike | None = None,
    head_count: int | None = None,
) -> Command:
    """Shuffle file."""
    cmd = Command("shuf")

    if head_count is not None:
        cmd += ["--head-count", head_count]

    if output is not None:
        cmd += ["--output", output]

    cmd += [filename]

    return cmd


def make_command_picard_markduplicate(
    *,
    aligned_read: PathLike,
    duped_read: PathLike,
    metrics: PathLike,
) -> Command:
    """Make command for mark duplicates.

    Args:
        aligned_read: SAM/BAM/CRAM file to analyze
        deduped_read: duplicated reads marked aligned reads
        metrics: duplicataion metrics

    Returns:
        command for picard markduplicates

    """
    return Command(
        "picard",
        "MarkDuplicates",
        "--INPUT",
        aligned_read,
        "--METRICS_FILE",
        metrics,
        "--OUTPUT",
        duped_read,
    )


def make_command_picard_collecthsmetrics(
    *,
    reference: PathLike,
    reference_dict: PathLike,
    aligned_read: PathLike,
    metrics: PathLike,
    bait_bed: PathLike,
    target_bed: PathLike,
    tmp_dir: PathLike,
    clip_overlapping_reads: bool = False,
    coverage_cap: int = 100000,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
    arguments_file: PathLike | None = None,
    base_coverage: PathLike | None = None,
    target_coverage: PathLike | None = None,
) -> Command:
    """Make command for picard CollectHsMetrics.

    Args:
        reference: reference genome in fasta
        reference_dict: reference sequence dictionary, where header contains sequence records
        aligned_read: aligned reead file in BAM
        metrics: metrics for hybrid-selection
        bait_bed: location of the baits in BED (to be converted into intervals)
        target_bed: location of the targets in BED (to be converted into intervals)
        tmp_dir: store intervals files
        clip_overlapping_reads: if set True, clip overlapping reads
        coverage_cap: max coverage limit for theoretical sensitivity calculation
        min_base_quality: min base quality to contribute to coverage
        min_mapping_quality: min mapping quality to contribute to coverage
        arguments_file: file with arguments, to be added to command line
        base_coverage: write per-base coverage information
        target_coverrage: write per-target coverage information

    Returns:
        command to run picard CollectHsMetrics

    """
    tmp_dir = parse_path(tmp_dir)
    bait = parse_path(bait_bed)
    target = parse_path(target_bed)

    bait_intervals = tmp_dir / bait.with_suffix(".intervals").name
    target_intervals = tmp_dir / target.with_suffix(".intervals").name

    bait_to_interval_cmd = make_command_picard_bedtointerval(bait, bait_intervals, reference_dict)
    target_to_interval_cmd = make_command_picard_bedtointerval(target, target_intervals, reference_dict)

    metrics_cmd = Command(
        "picard",
        "CollectHsMetrics",
        "--INPUT",
        aligned_read,
        "--OUTPUT",
        metrics,
        "--REFERENCE_SEQUENCE",
        reference,
        "--COVERAGE_CAP",
        coverage_cap,
        "--MINIMUM_BASE_QUALITY",
        min_base_quality,
        "--MINIMUM_MAPPING_QUALITY",
        min_mapping_quality,
        "--BAIT_INTERVALS",
        bait_intervals,
        "--TARGET_INTERVALS",
        target_intervals,
    )

    if clip_overlapping_reads:
        metrics_cmd += ["--CLIP_OVERLAPPING_READS"]

    if arguments_file is not None:
        metrics_cmd += ["--argument_file", arguments_file]

    if base_coverage is not None:
        metrics_cmd += ["--PER_BASE_COVERAGE", base_coverage]

    if target_coverage is not None:
        metrics_cmd += ["--PER_TARGET_COVERAGE", target_coverage]

    return bait_to_interval_cmd & target_to_interval_cmd & metrics_cmd


def make_command_picard_bedtointerval(
    bedfile: PathLike,
    intervals: PathLike,
    reference_dict: PathLike,
) -> Command:
    """Make command tp run picard BedToIntervalList.

    Args:
        bedfile: chr, start, stop coordinates in BED format
        intervals: coordinates converted into interval list
        reference_dict: reference sequence dictionary, where header contains sequence records

    Returns:
        command for picard bed to interval conversion.

    """
    return Command("picard", "BedToIntervalList", "-I", bedfile, "-O", intervals, "-SD", reference_dict)


def make_command_gatk_recalibrate(
    *,
    reference: PathLike,
    aligned_read: PathLike,
    calibrated_read: PathLike,
    calibrated_table: PathLike,
    known_sites: PathLike | SomePathLikes,
    intervals: PathLike | SomePathLikes | None = None,
) -> Command:
    """Make command for base quality score recalibrated aligned reads.

    Args:
        reference: reference genome in fasta
        aligned_read: aligned reads
        calibrated_read: recalibrated reads to used for variant calling
        calibrated_table: recalibration table file
        known_sites: one or more filepaths to known polymorphic sites to exclude
        intervals: one or more genomic intervals over which to operate

    Returns:
        command for GATK base quality score recalibration

    """
    recalibrate = Command(
        "gatk",
        "BaseRecalibrator",
        "-R",
        reference,
        "-I",
        aligned_read,
        "-O",
        calibrated_table,
    )
    for path in parse_paths(known_sites):
        recalibrate += ["--known-sites", path]

    if intervals is not None:
        for path in parse_paths(intervals):
            recalibrate += ["--intervals", path]

    apply = Command(
        "gatk",
        "ApplyBQSR",
        "-R",
        reference,
        "-I",
        aligned_read,
        "-O",
        calibrated_read,
        "--bqsr-recal-file",
        calibrated_table,
    )
    return recalibrate & apply


def make_command_gatk_haplotypecaller(
    *,
    reference: PathLike,
    aligned_read: PathLike,
    calls: PathLike,
    intervals: PathLike | SomePathLikes | None = None,
    no_soft_clipped_bases: bool = False,
    max_reads_per_start: int = 0,
    min_base_quality_score: int = 10,
    annotation_group: str | list[str] = "StandardAnnotation",
    padding: int | None = None,
    bamout: PathLike | None = None,
) -> Command:
    """Make command for GATK haplotypecaller.

    Args:
        reference: reference genome in fasta
        aligned_read: aligned reads
        calls: haplotypecaller calls
        intervals: one or more genomic intervals over which to operate
        no_soft_clipped_bases: if set, do not analyze soft clipped bases in the reads
        max_reads_per_start: maximum # of reads to retain per alignment start position. By default, set to 0 to disable.
        min_base_quality_score: minimum base quality required to consider a base for calling
        annotation_group: one or more groups of annotations to apply to variant calls
        padding: number of additional bases to include around each assembly region, if not set, 100 bp is used
        bamout: assembled haplotypes

    Returns:
        gatk haplotypecaller command

    """
    cmd = Command(
        "gatk",
        "HaplotypeCaller",
        "-R",
        reference,
        "-I",
        aligned_read,
        "-O",
        calls,
        "--min-base-quality-score",
        min_base_quality_score,
    )

    if intervals is not None:
        for path in parse_paths(intervals):
            cmd += ["--intervals", path]

    if no_soft_clipped_bases:
        cmd += ["--dont-use-soft-clipped-bases"]

    if max_reads_per_start is not None:
        cmd += ["--max-reads-per-alignment-start", max_reads_per_start]

    if bamout is not None:
        cmd += ["--bamout", bamout]

    annotation_group = [annotation_group] if isinstance(annotation_group, str) else annotation_group
    for group in annotation_group:
        cmd += ["-G", group]

    return cmd


def make_command_aggregate(
    *,
    script: PathLike,
    output_dir: PathLike,
    samplesheet: PathLike,
    bulk_db: PathLike | None = None,
    genes: PathLike | None = None,
    workflow: str,
) -> Command:
    """Make command for aggregate and summarize R script.

    # FIXME: aggregate and summarize portion of the R script --> python

    Args:
        script: path to script
        workflow: name of the workflow
        output_dir: output directory
        samplesheet: samplesheet
        snps_db: snp database
        genes: genecode gene listed in bed format

    Returns:
        command to execute R script

    """
    if workflow == "rna" and genes is None:
        raise ValueError(f"{workflow=} requires gene BED file")

    if workflow in ("sid", "tga") and bulk_db is None:
        raise ValueError(f"{workflow=} requires SNP database")

    cmd = Command(
        "Rscript",
        script,
        "--sheet",
        samplesheet,
        "--wdir",
        output_dir,
        "--workflow",
        workflow,
    )

    if bulk_db is not None:
        cmd += ["--bulk_db", bulk_db]

    if genes is not None:
        cmd += ["--genes", genes]

    return cmd


def make_command_rtg_vcfeval(
    *,
    baseline: PathLike,
    called_variants: PathLike,
    reference: PathLike,
    outdir: PathLike,
    regions: PathLike | None = None,
    evaluation_regions: PathLike | None = None,
    baseline_sample: str = "",
    calls_sample: str = "",
    roc_field: str = "",
    threads: int = 1,
) -> Command:
    """Make command for RTG tools vcfeval.

    FIXME: Enable to use explicit output filenames, not just top level output directory.
    RTG Tools vcfeval_ reports on comparisons between called variants and a baseline variant set.

    Args:
        baseline: VCF file containing baseline variants
        called_variants: VCF file containing called variants
        reference: reference in SDF-formatted fasta
        outdir: directory for output
        regions: if set, filepath to BED file of regions to read VCF records from that overlap
        evaluation_regions: if set, filepath to BED file of regions to evaluate
        baseline_sample: if set, sample name to use as baseline
        calls_sample: if set, sample name to use for calls
        roc_field: if set, alternative name to 'GQ' to use for the ROC score VCF format field
        threads: number of threads

    Returns:
        rtg vcfeval command

    Raises:
        ValueError: if baseline_sample is specified while calls_sample is not
        ValueError: if calls_sample is specified while baseline_sample is not

    .. _vcfeval:
        https://cdn.rawgit.com/RealTimeGenomics/rtg-tools/master/installer/resources/tools/RTGOperationsManual/rtg_command_reference.html#vcfeval

    """
    cmd = Command(
        "rtg",
        "vcfeval",
        "-b",
        baseline,
        "-c",
        called_variants,
        "-t",
        reference,
        "-o",
        outdir,
        "-T",
        threads,
    )

    if regions:
        cmd += [f"--bed-regions={regions}"]

    if evaluation_regions:
        cmd += ["-e", evaluation_regions]

    if baseline_sample or calls_sample:
        if not baseline_sample:
            raise ValueError("baseline sample name must be specified if test sample name is specified")
        elif not calls_sample:
            raise ValueError("test sample name must be specified if baseline sample name is specified")
        cmd += [f"--sample={baseline_sample},{calls_sample}"]

    if roc_field:
        cmd += ["-f", roc_field]

    return cmd


def make_command_cellranger_mkfastq(
    *,
    runid: str,
    run_dir: PathLike,
    ss: PathLike,
    output_dir: PathLike,
    delete_undetermined: bool = True,
    threads: int = 1,
) -> Command:
    """Make command for cellranger mkfastq.

    Args:
        runid: name of the folder created by mkfastq
        run_dir: Illumina run folder that contains the BCL files
        ss: samplesheet in comma-separated value (CSV) file
        output_dir: folder with FASTQs, reports and stats
        delete_undetermined: delete the Undetermined FASTQ file
        threads: max cores the pipeline may request at one time

    Returns:
        cellranger mkfastq command.

    """
    cmd = Command(
        "cellranger",
        "mkfastq",
        f"--id={runid}",
        f"--run={run_dir}",
        f"--csv={ss}",
        f"--output-dir={output_dir}",
        f"--localcores={threads}",
    )

    if delete_undetermined:
        cmd += ["--delete-undetermined"]

    return cmd


def make_command_cellranger_count(
    *,
    id: str,
    sample: str,
    transcriptome: PathLike,
    fastqs: PathLike,
    chemistry: str = "auto",
    use_intron: bool = True,
    threads: int = 1,
) -> Command:
    """Make command for cellranger count.

    Args:
        id: a unique run id and output folder name `[a-zA-Z0-9_-]+`
        sample: prefix of the filenames of FASTQs to select
        fastqs: path to input FASTQ data
        transcriptome: path of folder containing 10x-compatible transcriptome reference
        use_intron: if set True, include intronic reads in count. Default is True
        threads: set max cores the pipeline may request at one time. Only applies to local jobs
        chemistry: `auto` for autodetection, `SC3Pv1`, `SC3Pv2`, `SC3Pv3` for Single Cell 3' v1/v2/v3
    Returns:
        cellranger count command.

    """
    cmd = Command(
        "cellranger",
        "count",
        "--id",
        id,
        "--chemistry",
        chemistry,
        "--sample",
        sample,
        "--transcriptome",
        transcriptome,
        "--fastqs",
        fastqs,
        "--localcores",
        threads,
    )

    if not use_intron:
        cmd += ["--include-introns", "false"]

    return cmd


def make_command_cellranger_aggregate() -> Command:
    """Make command for cellranger aggr.

    Args:
        FIXME

    Returns:
        cellranger aggr command.

    """
    return Command()


def make_command_samtools_stats(
    *,
    aligned_read: PathLike,
    output_stats: PathLike,
    threads: int = 1,
) -> Command:
    """Make command for samtools stats.

    Args:
        aligned_read: input aligned read in BAM,SAM, CRAM
        output_read: output file name
        threads: number of threads

    Returns:
        samtools stats command

    """
    aligned_read = parse_path(aligned_read)
    output_stats = parse_path(output_stats)
    cmd = Command("samtools", "stats", "--threads", threads, aligned_read)

    return cmd.redirect_output(stdout=output_stats)


def make_command_read_counter(
    *,
    aligned_read: PathLike,
    output: PathLike,
    chromosome: str | None = None,
    window: int = 1000000,
    quality: int = 20,
) -> Command:
    """Make command for readCounter.

    Args:
        aligned_read: aligned read to be counted
        output: readCounter output wig file
        chromosome: chromosomes to be counted
        window: window size
        quality: quality threshold
        output: output wig file

    Returns:
        readCounter command from HMMCopy

    """
    if chromosome is None:
        chromosome = (
            "chr1,chr2,chr3,chr4,chr5,chr6,chr7,"
            "chr8,chr9,chr10,chr11,chr12,chr13,"
            "chr14,chr15,chr16,chr17,chr18,chr19,"
            "chr20,chr21,chr22,chrX,chrY"
        )

    if output is None:
        aligned_read = parse_path(aligned_read)
        output = aligned_read.parent / f"{aligned_read.stem}.wig"

    cmd = Command(
        "readCounter",
        "--window",
        window,
        "--quality",
        quality,
        "--chromosome",
        chromosome,
        aligned_read,
    )

    return cmd.redirect_output(stdout=output)


def make_command_run_ichorcna(
    *,
    sampleid: str,
    bigwig: PathLike,
    output: PathLike,
    ploidy: str | None = None,
    normal_fraction: str | None = None,
    include_homd: bool = False,
    chrs: str | None = None,
    chr_train: str | None = None,
    estimate_normal: bool = False,
    estimate_ploidy: bool = False,
    estimate_scprevalence: bool = False,
    txne: float = 0.9999,
    txn_strength: int = 10000,
    logfile: PathLike | None = None,
) -> Command:
    """Make command for runIchorCNA.R script.

    Args:
        sampleid: Sample id
        bigwig: BigWig file generated from aligned reads
        output: Output directory
        ploidy: Initial tumour ploidy
        normal_fraction: Initial normal contamination
        include_homd: If FALSE, then exclude HOMD state
        chrs: Specify chromosomes to analyze
        chr_train: Specify chromosomes to estimate params
        estimate_normal: Estimate normal
        estimate_ploidy: Estimate tumour ploidy
        estimate_scprevalence: Estimate subclonal prevalence
        txne: Self-transition probability
        txn_strength: Transition pseudo-counts
        logfile: log of ichorCNA

    Returns:
        readCounter command from HMMCopy

    """
    if ploidy is None:
        ploidy = "c(1.75,2,2.25)"

    if normal_fraction is None:
        normal_fraction = "c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)"

    if chrs is None:
        chrs = 'c(1:22,"X")'

    if chr_train is None:
        chr_train = "c(1:22)"

    cmd = Command(
        "Rscript",
        parse_path("./askcell/dna/runIchorCNA.R"),
        "--id",
        sampleid,
        "--WIG",
        parse_path(bigwig),
        "--ploidy",
        ploidy,
        "--normal",
        normal_fraction,
        "--gcWig",
        parse_path("./askcell/data/gc_hg38_1000kb.wig"),
        "--mapWig",
        parse_path("./askcell/data/map_hg38_1000kb.wig"),
        "--centromere",
        parse_path("./askcell/data/GRCh38.GCA_000001405.2_centromere_acen.txt"),
        "--normalPanel",
        parse_path("./askcell/data/cnv.pon_median.rds"),
        "--chrs",
        chrs,
        "--chrTrain",
        chr_train,
        "--txnE",
        txne,
        "--txnStrength",
        txn_strength,
        "--outDir",
        parse_path(output),
    )

    if include_homd:
        cmd += ["--includeHOMD", True]

    if estimate_normal:
        cmd += ["--estimateNormal", True]

    if estimate_ploidy:
        cmd += ["--estimatePloidy", True]

    if estimate_scprevalence:
        cmd += ["--estimateScPrevalence", True]

    return cmd.redirect_output(stdout=logfile)
