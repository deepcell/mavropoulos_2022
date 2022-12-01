==============================================
Add demux metrics as part of run-level metrics
==============================================

Description:

This is 1 of 3 PRs addressing QC metrics for askcell workflow.

Input is sequencing run directory (Illumina generated outputs) and collect the following run-level metrics.

Output is run-level results in DataFrame format (2 rows, where 1st row is a header and 2nd row contains the value in int or float format).

* Total yield: this is total yield from a run

* >Q30%: percentage of bases with a quality score of 30 or higher

* % of Clusters PF: percentage of clusters passing Illumina Chastity filter (clear signal in the base)

* Density of clusters on the flowcell in thousands per mm^2 (not used for Novaseq, HiSeq 4000 and X) to evaluate over or under-loading of the library

Stories:

    1. `BIOIN-275 <https://deepcellbio.atlassian.net/browse/BIOIN-275>`_: Add demux metrics as part of run-level metrics
