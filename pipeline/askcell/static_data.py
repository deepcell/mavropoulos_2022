"""static_data.py static_data to go into configuration.

NOTE: No need to use ALL CAPS naming convention here since all contents of this module are global data.

TODO: Move to configuration

"""
import os

from .utils.pathutils import parse_path


# bin
PICARD = parse_path(os.environ.get("PICARD", "/opt/miniconda3/envs/workflow/share/picard-2.26.2-0/picard.jar"))

# R script
# FIXME: Remove all uses of __file__ ASAP!!!!
R_SUMMARIZE = parse_path(__file__).parent / "summarize.R"
R_WGS_SUMMARIZE = parse_path(__file__).parent / "dna" / "wgs_postprocessing_titan2.R"

# Human
BASE = parse_path(os.environ.get("BASE", "/mnt/ro/BioPipeline/resource/Hsapien/GRCh38"))
GENOME = BASE / "genome/GRCh38_no_alt.p13.genome.fa"
GENOME_DICT = BASE / "genome/GRCh38_no_alt.p13.genome.dict"
GENOME_SDF = BASE / "genome/sdf_no_alt"
rRNA_FA = BASE / "index/bwa_rRNA/rRNA.fasta"
GENCODE_GENE_BED = BASE / "annotation/gencode.v36.annotation.genes.bed6"
GENOME_INDEX = BASE / "index/bwa/GRCh38_no_alt.p13.genome.fa"
STAR_GENOME_INDEX = BASE / "index/star_151bp"
GNOMAD = BASE / "gnomad_v3/af-only-gnomad.hg38.vcf.gz"
WGS_BASELINE = parse_path(
    os.environ.get("WGS_BASELINE", "/mnt/ro/Biology/NGS/WGS/WGS14/analysis_files/genome_mb_coverage.csv")
)

# Mouse
MM_BASE = parse_path(os.environ.get("MM_BASE", "/mnt/ro/BioPipeline/resource/Mouse/GRCm39"))
MM_rRNA_FA = MM_BASE / "index/bwa_rRNA/mm_rRNA.fasta"
STAR_MM_GENOME_INDEX = MM_BASE / "index/mm_star_151bp"

# SWIFT platform-specific
PANELS = BASE / "panels"
SWIFT_SID_SNP = PANELS / "swift_sid/swift_sid.snp.grch38.bed"
SWIFT_SID_AMPLICON = PANELS / "swift_sid/swift_sid.amplicons.grch38.bed"
SWIFT_SID_PRIMER = PANELS / "swift_sid/swift_sid.primers.grch38.bed"
SWIFT_LUNG_AMPLICON = PANELS / "swift_lung/swift_lung.amplicons.grch38.bed"
SWIFT_LUNG_TARGET = PANELS / "swift_lung/swift_lung.targets.grch38.bed"
SWIFT_LUNG_PRIMER = PANELS / "swift_lung/swift_lung.primers.grch38.bed"
SWIFT_72G_AMPLICON = PANELS / "swift_72g/swift_72g.amplicons.grch38.bed"
SWIFT_72G_TARGET = PANELS / "swift_72g/swift_72g.targets.grch38.bed"
SWIFT_72G_PRIMER = PANELS / "swift_72g/swift_72g.primers.grch38.bed"
SWIFT_BRCA_AMPLICON = PANELS / "swift_brca/swift_brca.amplicons.grch38.bed"
SWIFT_BRCA_TARGET = PANELS / "swift_brca/swift_brca.targets.grch38.bed"
SWIFT_BRCA_PRIMER = PANELS / "swift_brca/swift_brca.primers.grch38.bed"
SWIFT_SID_SNP_DB = BASE / "pon/bulk_pure_cell_line_swift_sid/bulk_pure_cell_line_swift_sid.snp.rds"
AMPLISEQ_HOTSPOT2_TARGET = PANELS / "ampliseq_hotspot2/ampliseq.targets.hg38.bed"
AMPLISEQ_HOTSPOT2_PRIMER = PANELS / "ampliseq_hotspot2/ampliseq.primers.hg38.bed"
AMPLISEQ_HOTSPOT2_AMPLICON = PANELS / "ampliseq_hotspot2/ampliseq.amplicons.hg38.bed"
TWIST_EXOME_TARGET = PANELS / "twist_exome/hg38_Twist_exome_2_1_annotated_targets.bed"

# GIAB
GIAB = BASE / "goldset/giab/GIAB/"
NA12878_CALLS = GIAB / "NA12878_HG001/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
NA12878_CONFIDENT_REGIONS = GIAB / "NA12878_HG001/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"

# Universal primer
P5 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
P7 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"

# Panel to target BED mapping
PANELS_TARGET_KEY = {
    "swift_sid": SWIFT_SID_SNP,
    "swift_lung": SWIFT_LUNG_TARGET,
    "swift_brca": SWIFT_BRCA_TARGET,
    "swift_72g": SWIFT_72G_TARGET,
    "ampliseq_hotspot2": AMPLISEQ_HOTSPOT2_TARGET,
}

# FIXME: collect_metrics use this file, not Amplicon files
PANELS_AMPLICON_KEY = {
    "swift_sid": SWIFT_SID_AMPLICON,
    "swift_lung": SWIFT_LUNG_TARGET,
    "swift_brca": SWIFT_BRCA_TARGET,
    "swift_72g": SWIFT_72G_TARGET,
    "ampliseq_hotspot2": AMPLISEQ_HOTSPOT2_TARGET,
}

PANELS_PRIMER_KEY = {
    "swift_sid": SWIFT_SID_PRIMER,
    "swift_lung": SWIFT_LUNG_PRIMER,
    "swift_brca": SWIFT_BRCA_PRIMER,
    "swift_72g": SWIFT_72G_PRIMER,
    "ampliseq_hotspot2": AMPLISEQ_HOTSPOT2_PRIMER,
}

# rnaseq trim parameters
AG_TRIM = [
    "A{5}G{5}",
    "A{10}G{10}",
    "A{5}G{10}",
    "A{10}G{5}",
    "A{15}G{15}",
    "A{15}",
    "A{20}",
    "A{30}",
    "A{40}",
    "A{50}",
    "A{60}",
    "G{20}",
    "G{30}",
    "G{40}",
]

REFERENCE_CONTROLS_CALLS = {"NA12878": NA12878_CALLS}
REFERENCE_CONTROLS_CONFIDENT_REGIONS = {"NA12878": NA12878_CONFIDENT_REGIONS}
