"""Stat Reader Module.

This module provides support for reading, parsing, and summarizing import metrics produced
during the demultiplexing of an Illumina sequencing run.

"""

from __future__ import annotations

import statistics
import xml.etree.ElementTree as ET

from collections import defaultdict
from typing import Any, Sequence

from . import PathLike
from .utils.fileutils import load_json


def read_conversionstat(filename: PathLike) -> dict[str, Any]:
    """Read ConversionStats.xml file Illumina output, retrieving some basic demultiplexing metrics.

    ConversionStats.xml contains information on basecalling per tile including
    raw and pf cluster counts.
    Output dictionary contains `Raw` and `Pf` values for each sample for each lane.

    Cluster count statistics include raw and passing filter cluster information.

    ======= ========== ======================================
    Name    Type       Definition
    ======= ========== ======================================
    Raw     int        count of raw clusters
    Pf      int        count of clusters passing filter
    ======= ========== ======================================

    Note:
        Input file is in Stats directory in Illumina sequencing output directory.

    Args:
        filename: full path to ConversionStats.xml file.

    Returns:
        dictionary of sample-specific raw and passing filter cluster count.

    """
    convstats: dict[str, Any] = defaultdict(dict)
    root = ET.parse(filename).getroot()

    flowcell = root.find("Flowcell")

    if not flowcell:
        raise ValueError("invalid conversion stats")

    convstats["FCID"] = flowcell.attrib["flowcell-id"]

    for sample in root.findall("Flowcell/Project/Sample"):
        sid = sample.attrib["name"]
        if sid == "all":
            continue
        for barcode in sample.findall("Barcode"):
            if barcode.attrib["name"] != "all":
                for lane in barcode.findall("Lane"):
                    lane_key = f'LANE_{lane.attrib["number"]}'
                    if lane_key not in convstats[sid]:
                        convstats[sid][lane_key] = {"Raw": 0, "PF": 0}
                    for raw in lane.findall("./Tile/Raw"):
                        convstats[sid][lane_key]["Raw"] += int(raw.findtext("ClusterCount") or 0)
                    for pf in lane.findall("./Tile/Pf"):
                        convstats[sid][lane_key]["PF"] += int(pf.findtext("ClusterCount") or 0)

    return convstats


def read_demuxstat(filename: PathLike, convstats: dict[str, Any]) -> dict[str, Any]:
    """Read Stats.json and call read_conversionstat and return dictionary of demultiplexing statistics.

    Stats.json is default output from demultiplexing containing overall summary statistics, found
    in Illumina sequencing output Stats directory.
    convstats is an output from :func:`parse_conversionstat`.

    Following metrics reported for each sample identifier.

    ============================ ======= =====================================================
    Name                         Type    Definition
    ============================ ======= =====================================================
    CONTROL                      str     control_type if available else null
    INDEX                        tuple   index in form of tuple of str
    PF_Q30_BASES_PCT             float   % bases with qscore > 30 from clusters passing filter
    INDEX_READS_ONE_MISMATCH_PCT float   % of indexes with 1 mismatch
    NUM_READS                    int     number of reads
    PF_PCT                       float   % of clusters passing filter
    YIELD_MBASES                 int     yielded bases in MB
    PF_QUALITY_SCORE_MEAN        float   quality score mean
    CLUSTERS_RAW                 int     number of clusters
    CLUSTERS_PF                  int     number of clusters passing filters
    CLUSTERS_RAW_PER_LANE_PCT    float   % raw reads in given sid over total
    INDEX_READS_PERFECT_PCT      float   % of indexes with 0 mismatch
    ============================ ======= =====================================================

    Args:
        filename:  path to Stats.json summary statistic output
        convstats: a dictionary containing raw and pf clusters count information

    Returns:
        information from Table 1 for each sample identifier

    """
    demuxstats_keys = [
        "CONTROL",
        "INDEX",
        "PF_Q30_BASES_PCT",
        "INDEX_READS_ONE_MISMATCH_PCT",
        "NUM_READS",
        "PF_PCT",
        "YIELD_MBASES",
        "PF_QUALITY_SCORE_MEAN",
        "CLUSTERS_RAW",
        "CLUSTERS_PF",
        "CLUSTERS_RAW_PER_LANE_PCT",
        "INDEX_READS_PERFECT_PCT",
    ]

    input_data = load_json(filename)

    demuxstats: dict[str, Any] = {}
    conversion_results = input_data["ConversionResults"]
    for lane in conversion_results:
        lane_key = "LANE_%s" % lane["LaneNumber"]
        demuxstats[lane_key] = {}
        demux_result = lane["DemuxResults"]
        demux_result.append(lane["Undetermined"])

        for result in demux_result:
            sid = result.get("SampleName", "Undetermined")
            sid_convstats = convstats[sid][lane_key]

            dstat: dict[str, Any] = {key: None for key in demuxstats_keys}
            dstat["CONTROL"] = None
            dstat["NUM_READS"] = result["NumberReads"]
            dstat["CLUSTERS_RAW"] = sid_convstats["Raw"]
            dstat["CLUSTERS_PF"] = sid_convstats["PF"]
            dstat["CLUSTERS_RAW_PER_LANE_PCT"] = (
                100.0 * sid_convstats["Raw"] / lane["TotalClustersRaw"] if (lane["TotalClustersRaw"]) else 0
            )
            dstat["PF_PCT"] = 100.0 * sid_convstats["PF"] / sid_convstats["Raw"] if (sid_convstats["Raw"]) else 0

            if sid != "Undetermined":
                indexmetrics = result["IndexMetrics"]
                mm0 = sum(read["MismatchCounts"]["0"] for read in indexmetrics)
                mm1 = sum(read["MismatchCounts"]["1"] for read in indexmetrics)
                dstat["INDEX"] = sorted(read["IndexSequence"] for read in indexmetrics)
                dstat["INDEX_READS_ONE_MISMATCH_PCT"] = 100.0 * mm1 / (mm0 + mm1) if (mm0 + mm1) else 0
                dstat["INDEX_READS_PERFECT_PCT"] = 100.0 * mm0 / (mm0 + mm1) if (mm0 + mm1) else 0

            readmetrics = result["ReadMetrics"]
            q30_base = sum(read["YieldQ30"] for read in readmetrics)
            total_base = sum(read["Yield"] for read in readmetrics)
            quality_score_sum = sum(read["QualityScoreSum"] for read in readmetrics)

            dstat["PF_Q30_BASES_PCT"] = 100.0 * q30_base / total_base if total_base else 0
            dstat["YIELD_MBASES"] = round(total_base / 1000000)
            dstat["PF_QUALITY_SCORE_MEAN"] = quality_score_sum / total_base if total_base else 0
            demuxstats[lane_key][sid] = dstat

    return demuxstats


def summarize_demuxstat(data: dict[str, Any]) -> dict[str, Any]:
    """Summarize gathered demultiplexing statistics.

    Sample-level data is agglomerated into a single run-level metrics.
    Run-level total data is summed across samples across all lanes.

    Args:
        data: gathered data during demultiplexing steps

    Returns:
        dictionary of agglomerated statistics of the form::

         {
          'CLUSTERS_RAW_TOTAL': int, total number of clusters,
          'CLUSTERS_PF_TOTAL': int, total number of clusters passed filter,
          'READS_TOTAL': int, total yielded indexed and unindexed reads,
          'READS_PER_SAMPLE_MEDIAN': float, median yielded read per indexed sample,
          'READS_UNINDEXED_TOTAL': int, total number of unindexed reads,
          'READS_UNINDEXED_PCT': float, percentage of unindexed reads over total reads
         }

    """
    composite: dict[str, Any] = {}
    composite["CLUSTERS_RAW_TOTAL"] = sum(sid["CLUSTERS_RAW"] for lane in data.values() for sid in lane.values())
    composite["CLUSTERS_PF_TOTAL"] = sum(sid["CLUSTERS_PF"] for lane in data.values() for sid in lane.values())
    composite["READS_TOTAL"] = sum(sid["NUM_READS"] for lane in data.values() for sid in lane.values())
    composite["READS_UNINDEXED_TOTAL"] = sum(lane["Undetermined"]["NUM_READS"] for lane in data.values())
    composite["READS_UNINDEXED_PCT"] = (
        100.0 * composite["READS_UNINDEXED_TOTAL"] / composite["READS_TOTAL"] if (composite["READS_TOTAL"]) else 0
    )

    reads: list[int] = []
    for lane in data:
        reads += get_list_from_samples(data[lane], "NUM_READS", sids=get_sids(data[lane], excludes=["Undetermined"]))
    composite["READS_PER_SAMPLE_MEDIAN"] = statistics.median(reads)

    return composite


def get_sids(data: dict[str, Any], excludes: Sequence[str] | None = None) -> list[str]:
    """Return list of sids.

    Args:
        data: dictionary containing sid as key.
        excludes: if not None, provide all keys in data dictionary. \
            sids in excludes list should be excluded.

    Returns:
        list of sids

    """
    return list(data) if not excludes else [sid for sid in data if sid not in excludes]


def get_list_from_samples(data: dict[str, Any], key: str, sids: Sequence[str] | None = None) -> list[int]:
    """Return list of values noted by key.

    Args:
        data: dictionary of values
        key:  key for dictionary
        sids: if not None, include all key values from all samples

    Returns:
        list of data

    """
    return [data[sid][key] for sid in (sids or data)]
