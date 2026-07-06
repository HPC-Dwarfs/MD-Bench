* Evaluate LJ mixing rules on GPUs and other cases
* Evaluate CP for AMD GPUs
* Evaluate Lennard-Jones (and Coloumb) force components to be integrated into short-range kernels
* SetupPBC should run on GPU (Verlet Lists this is clear, check Cluster Pair)
  - Measure its impact (apparently its contribution goes to the Neighbor part)
* Evaluate performance of SIMD compressed kernels and see if they can improve
* Implement LJ case from https://ieeexplore.ieee.org/document/11370954 for ARM and x86
* Implement stubbed case and gather benchmark for GPUs
