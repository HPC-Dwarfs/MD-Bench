* Evaluate LJ mixing rules on GPUs and other cases
* Evaluate CP for AMD GPUs
* Use a single capacity for the neighbor-lists and evaluate CPU vs GPU performance
* Evaluate Lennard-Jones (and Coloumb) force components to be integrated into short-range kernels
* Double cut-off method with pruning (inner, outer): controlled by the additive `outer_skin` parameter (outer cutoff = cutforce + skin + outer_skin; default 0 disables the scheme, no separate toggle). CPU clusterpair done (partition prune in src/clusterpair/neighbor.c; force kernels use inner counts in src/clusterpair/force_lj.c). GPU prune kernels run device-resident (cudaPruneNeighbor / cudaPruneNeighborSup in src/clusterpair/force_lj_cuda{,_sup}.cu), dispatched from pruneNeighborCUDA. Follow-ups:
  - Pad inner-section tail for CLUSTER_N < VECTOR_WIDTH (2xNN); v1 works but may visit extra j-clusters at the tail (filtered by per-pair cutforcesq, so still correct)
  - Extend to verletlist (no prune infrastructure today)
* SetupPBC should run on GPU (Verlet Lists this is clear, check Cluster Pair)
  - Measure its impact (apparently its contribution goes to the Neighbor part)
* Add displacement-based strategy to trigger neighbor-lists rebuilding
* Implement compression of atoms that need to be computed, only execute arithmetic when register is full
* Implement LJ case from https://ieeexplore.ieee.org/document/11370954 for ARM and x86
* Implement stubbed case and gather benchmark for GPUs
