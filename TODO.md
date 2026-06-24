* Evaluate LJ mixing rules on GPUs and other cases
* Evaluate CP for AMD GPUs
* Use a single capacity for the neighbor-lists and evaluate CPU vs GPU performance
* Evaluate Lennard-Jones (and Coloumb) force components to be integrated into short-range kernels
* SetupPBC should run on GPU (Verlet Lists this is clear, check Cluster Pair)
  - Measure its impact (apparently its contribution goes to the Neighbor part)
* Add displacement-based strategy to trigger neighbor-lists rebuilding
* Implement compression of atoms that need to be computed, only execute arithmetic when register is full
* Implement LJ case from https://ieeexplore.ieee.org/document/11370954 for ARM and x86
* Implement stubbed case and gather benchmark for GPUs
