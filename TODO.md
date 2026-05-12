* Apply Lennard-Jones Geometric Combination Rule
* Define frequency to sort atoms (independent neighboring frequency)
* Use a single capacity for the neighbor-lists and evaluate CPU vs GPU performance
* Evaluate Lennard-Jones (and Coloumb) force components to be integrated into short-range kernels
* Double cut-off method with pruning (inner, outer)
* Implement compression of atoms that need to be computed, only execute
arithmetic when register is full
* Implement LJ case from https://ieeexplore.ieee.org/document/11370954 for ARM and x86
* Implement stubbed case and gather benchmark for GPUs
