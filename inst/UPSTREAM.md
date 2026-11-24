# SpatialPCA attribution

Refinement adapted from `refine_cluster_10x()` in
https://github.com/shangll123/SpatialPCA/blob/1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66/R/SpatialPCA_utilties.R

Commit: `1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66`.
Original header: Package SpatialPCA; version 1.1.0; date 2021-10-27;
authors L. Shang and X. Zhou; University of Michigan, Department of Biostatistics.
Repository DESCRIPTION at that commit declares GPL (>= 3). Copies of the source
and DESCRIPTION are retained in `inst/upstream/` with their original notices.
MSGPCA is distributed under GPL (>= 3); see the included GPL-3 text.

Changes: renamed function; character label-name lookup fixes numeric/factor
indexing; explicit validation; rejection of duplicate coordinates and insufficient
sample sizes; exact k-NN/radius queries retrieve boundary ties, followed by exact Euclidean
distances and input-order sorting. A rowwise reference engine remains available.
Both avoid the full distance matrix.
The single sequential pass, focal inclusion, thresholds, and ordering ties by
observation index are preserved. This is not an iterative majority-vote smoother.
No SpatialPCA package dependency is introduced.

Citation: Shang, L. and Zhou, X. (2022). Spatially aware dimension reduction for
spatial transcriptomics. Nature Communications 13, 7203.
https://doi.org/10.1038/s41467-022-34879-1

Code redistribution here does not imply upstream endorsement. Authorship and
public release authority for MSGPCA remain with the manuscript's investigators.
