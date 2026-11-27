* Evaluate Lennard-Jones (and Coloumb) force components to be integrated into short-range kernels
* SetupPBC should run on GPU (Verlet Lists this is clear, check Cluster Pair)
  - Measure its impact (apparently its contribution goes to the Neighbor part)
* Evaluate performance of SIMD compressed kernels and see if they can be improved
* Implement LJ case from https://ieeexplore.ieee.org/document/11370954 for ARM and x86
* Evaluate stubbed case for GPUs and implement super-cluster case
* Use the same "GPU slicing" idea with MPI to improve strong scaling
  - Intra-node (OpenMP) j-slicing: split long cluster neighbor lists into (ci, j_begin, j_end) work items,
    parallelize over them and reduce i-forces with per-thread buffers or atomics (no comm.c changes)
  - Boundary pair ownership balancing (halfShell/eightShell): dynamically choose which of the two ranks
    computes a cross-boundary pair (both already hold the data, reverse comm returns the forces), driven by
    the balance.h timers; only helps inhomogeneous systems
