# MD-Bench

| :exclamation:  The repository was restructured and cleaned up and in this process the history was rewritten. Please make fresh forks and clones.   |
|-----------------------------------------|

MD-Bench is a toolbox for the performance engineering of short-range force
calculation kernels on molecular-dynamics applications. It aims at covering
state-of-the-art algorithms from different community codes such as
LAMMPS and GROMACS.

## Getting started

Clone the repository from GitHub:

```shell=
git clone https://github.com/RRZE-HPC/MD-Bench.git
```

Edit config.mk and configure the compiler toolchain

```makefile=
# Compiler tool chain (GCC/CLANG/ICC/ICX/ONEAPI/NVCC/HIPCC)
TOOLCHAIN ?= GCC
```

Best supported are ICC (deprecated legacy Intel compiler) and ICX (LLVM based
Intel compiler). Choose NVCC to enable CUDA GPU kernels (Nvidia), or HIPCC to
enable HIP GPU kernels (AMD).

The toolchain settings are located in the `./make` directory. Review the
settings for the configured toolchain. You can configure different settings in
`config.mk`, for starters on a X86 based system the defaults are fine.

To build the binary call, don't forget to load the compiler module on the
NHR@FAU clusters (e.g. `module load intel`):

```shell=
make
```

While the Makefile works with any version of GNU make, some features require GNU
make > v4.

You can run MD-Bench without any arguments:

```shell=
./MDBench-VL-ICX-X86-AVX2-DP
```

## Build system

MD-Bench uses a Makefile with pattern rules and automatic dependency generation.
If you add source files you do not need to change the Makefile as long as the
sources are placed either in the `./src/verletlist/`, `./src/clusterpair/` or
`./src/common` directories. If you change a file, all object files that depend
on it are rebuild.

All configuration variables can be overwritten from the command line, e.g. to
build with ICC without changing `./config.mk` build with:

```shell=
make TOOLCHAIN=ICC
```

Multiple configurations can be build at the same time. Every configuration has a
unique binary name `./MDBENCH-<build tag>`. Intermediate build results are
located in a `./build/build-<build tag>/` directory.

All make targets act on the current configuration set in `./config.mk`, but this
can be of course overwritten on the command line.

Supported make targets:

- `make`: Build the binary for current configuration.
- `make clean`: Remove intermediate build results.
- `make distclean`: Remove intermediate build results and binary. Also removes
generated tags and clangd files, more on that later.
- `make cleanall`: Remove all generated files. **Note**: This target applies to
all configurations.
- `make info` Output compiler version, useful for logging in automated benchmark
scripts.
- `make asm`: Generate assembly output of all source files. The assembly files
are placed in the intermediate build directory.
- `make format`: Reformat all source files with `clang-format` using the format
specification in `.clang-format`.

### Build time options

- `TOOLCHAIN`: Determines which toolchain makefile is included
(`GCC`/`CLANG`/`ICC`/`ICX`/`ONEAPI`/`NVCC`/`HIPCC`)
- `ISA`: No usage apart from tag strings
- `SIMD`: Controls the generation of intrinsic kernels for clusterpair.
X86 options: `NONE`/`SSE`/`AVX`/`AVX_FMA`/`AVX2`/`AVX512`.
ARM options: `NONE`/`NEON`/`SVE`/`SVE2`.
- `OPT_SCHEME`: Algorithmic variant (`verletlist` or `clusterpair`), different
source directories and main routines are used
- `ENABLE_LIKWID`: Turn on LIKWID instrumentation, the LIKWID library has to be available
- `ENABLE_OPENMP`: Enable OpenMP thread-level parallelization (true or false)
- `ENABLE_MPI`: Turn on the MPI parallel version of the code
- `DATA_TYPE`: Switch between single precision (`SP`) and double precision (`DP`) floating
point. This is controlled by defines.
- `ATOM_DATA_LAYOUT`: Switch between array-of-structure (`AOS`) and structure-of-array
(`SOA`) layout for atom positions and forces. Tradeoff between better cache
utilisation and easier SIMD vectorization.
- `NBLIST_DATA_LAYOUT`: Neighbor-list data layout (`auto`/`AOS`/`SOA`). `auto` defaults
to `AOS` for CPU builds and `SOA` for GPU builds.
- `LJ_COMB_RULE`: Lennard-Jones combination rule. `single` uses a single atom type
with broadcast global parameters (fastest, no type lookup). `geometric` uses
per-type parameters with geometric combination (default). `full` uses a full
type-pair matrix lookup (not supported in SIMD kernels).
- `COULOMB_RULE`: Electrostatics for the SIMD cluster-pair kernels (4xN, 2xNN, 2xN;
`none` (default)/`rf`/`ewald`). `rf` is a reaction-field with infinite dielectric, `ewald`
the real-space Ewald part with the analytical force correction used by GROMACS. The
reference and GPU kernels, the verletlist scheme and MPI runs do not support it.
Atom charges are given per atom type.
- `EXCLUSIONS` (`auto` (default)/`true`/`false`): Molecular exclusions. Atoms of the
same molecule (consecutive atom ids, `atom->molsize` atoms per molecule) do not
interact via LJ or the direct Coulomb term, but get the exclusion correction of the
selected `COULOMB_RULE`. `auto` enables it whenever `COULOMB_RULE` is not `none`.
- `ACCURATE_RECIPROCAL` (default: false): Refine the hardware reciprocal and
reciprocal square root estimates (12-14 bits) used by the cluster-pair LJ kernels with
Newton-Raphson steps. Needed
when forces have to agree with a code that evaluates 1/r^2 to full precision, e.g. to
compare against GROMACS.
- `DEBUG`: Enable additional debug output
- `SORT_ATOMS`: Resort atoms to ensure that atoms that are nearby are also close
to each other in the data structures
- `USE_REFERENCE_KERNEL`: Enforce usage of the C reference implementation for
validation (automatically set for GPU and `SIMD=NONE` builds)
- `USE_SCALAR_KERNEL`: Use scalar version and rely on the compiler for
auto-vectorization
- `USE_SIMD_KERNEL`: Use hand-written SIMD intrinsic kernels for force computation
- `USE_CUDA_HOST_MEMORY`: Enable pinned host memory for faster host-device transfers
- `USE_PINNED_MEMORY` (default: true): Allocate the host-side buffers that are
transferred to/from the GPU every reneighbor step (positions/velocities,
per-type LJ params, neighbor lists) as pinned memory instead of regular
pageable memory, for faster/more predictable `cudaMemcpy` transfers
- `GPU_ARCH`: Target GPU architecture for HIPCC (AMD) builds (e.g. `gfx90a`).
Must be set explicitly; common values: `gfx90a` (MI250X), `gfx940` (MI300A),
`gfx942` (MI300X).
- `CLUSTER_PAIR_KERNEL`: Kernel variant for the clusterpair scheme
(`auto`/`4xN`/`2xN`/`2xNN`/`gpusimple`/`supercluster`)
- `SUPERCLUSTER_DATA_LAYOUT`: Data layout for super-clustering kernels
(`AOS3`/`AOS4`/`SOA`)
- `SUPERCLUSTER_INVERSE_THREAD_MAPPING`: Map `threadIdx.y` to `cii` and
`threadIdx.x` to `cjj` (true or false). If false, uses the same thread
mapping and reduction as GROMACS.
- `SUPERCLUSTER_J_SLICES` (default: 4): Number of thread-block slices the GPU
grid is split into along each i-supercluster's j-neighbor list, for the
super-clustering kernel. 1 means each i-supercluster's neighbors are handled
by a single block; >1 splits the neighbor list across multiple concurrent
blocks (more parallelism/occupancy) that combine partial forces with
atomicAdd (more atomic contention). Tune based on system size and occupancy.
- `XTC_OUTPUT`: Enable XTC trajectory output (GROMACS file format)
- `MEM_TRACER`: Trace memory addresses for cache simulation (true or false)
- `INDEX_TRACER`: Trace indexes and distances for gather-MD analysis (true or false)
- `COMPUTE_STATS`: Compute and print additional runtime statistics

### Build for GPU targets

MD-Bench supports GPU kernels for both Nvidia (CUDA) and AMD (HIP) GPUs.
A binary either supports CPU kernels or GPU kernels — not both at the same time.

**Nvidia (CUDA):** Specify `NVCC` as toolchain. CUDA source files share the same
directories as CPU sources, using a `Cuda` suffix and `.cu` file extension. All
supported kernels are automatically switched to their CUDA variants at build time.

```shell=
make TOOLCHAIN=NVCC
```

**AMD (HIP/ROCm):** Specify `HIPCC` as toolchain and set `GPU_ARCH` to your
target architecture. `GPU_ARCH` has no default and must always be provided:

```shell=
make TOOLCHAIN=HIPCC GPU_ARCH=gfx90a
```

Common `GPU_ARCH` values: `gfx90a` (MI250X), `gfx940` (MI300A APU),
`gfx942` (MI300X).

## Command line arguments

MD-Bench can be executed without any arguments, in this case the full neighbor
list testcase with LJ force will be computed for 200 steps and a size of
32x32x32 unit cells.

- `-p / --params <string>`: file to read parameters from (can be specified more
than once). Default initialization sets parameters for default LJ testcase.
Supports `outer_skin = <value>` to enable double-cutoff pruning (outer cutoff =
`cutforce + skin + outer_skin`; default 0.0 disables it).
- `-f <string>`: force field (`lj`, `eam`, `lj_table`), default `lj`. For
anything different than `lj` you also need to provide a specific parameter file.
`lj_table` uses tabulated cubic Hermite spline interpolation. It is available for
`OPT_SCHEME=verletlist`, and for `OPT_SCHEME=clusterpair` only when built with the
reference kernel (`USE_REFERENCE_KERNEL=true`).
- `-i <string>`:  input file with atom positions (dump). MD-Bench supports
Brookhaven protein data bank (.pdb), GROMACS GROMOS87 (.gro), and LAMMPS dump
(.dmp) file formats
- `-e <string>`:  input file for EAM parameters
- `-n / --nsteps <int>`:  set number of timesteps for simulation (default 200)
- `-nx/-ny/-nz <int>`:  set linear dimension of systembox in x/y/z direction
(default 32 in every dimension)
- `-half <int>`: use half (1) or full (0) neighbor lists (default 0 - full
neighbor list). The cluster-pair SIMD kernels without MPI keep only half of
the periodic ghost clusters with half lists, so a half list has about half the entries of
a full list. Half lists are not thread safe, use one OpenMP thread
- `-r / --radius <real>`:   set cutoff radius (default 2.5)
- `-s / --skin <real>`:   set skin (verlet buffer, default 0.3)
- `-w <file>`:  write input atoms to file
- `--freq <real>`:  processor frequency (GHz), used to calculate cycle metrics
(default 2.4)
- `--vtk <string>`:    VTK output file for visualization
- `-gmxbench`:  use the water system of the GROMACS non-bonded benchmark
(`src/gromacs/nbnxm/benchmark`) instead of the default testcase: 3000 SPC/E-like
atoms with charges and intramolecular exclusions, cut-off 1.0 and no skin unless
`-r`/`-s` are given. `-size <int>` replicates the box (power of two) as `gmx
nonbonded-benchmark -size` does. Build with `COULOMB_RULE=ewald` (or `rf`),
`CLUSTER_PAIR_KERNEL=4xN` and `ACCURATE_RECIPROCAL=true` to compare with GROMACS.
- `--kernel-bench <int>`:  after setup only run the force kernel this many times on the
fixed neighbor lists and report the time per call and the useful pair rate, as the
GROMACS benchmark does
- `--dump-forces <file>`:  write the forces of a single kernel pass (`atom id fx fy
fz`, original atom order) to a text file

## Available testcases

For all variants you can switch between single precision and double precision
and between AOS versus SOA data layouts using build time options. You can use
the half neighbour list algorithm instead of the default full neighbour list by
setting `-half 1`.

### Lennard-Jones potential for solid copper

Just start without any command line argument, this is the default testcase. You
may change the number of timesteps using the `-n` options and change the problem
size using the `-nz, -ny, -nz` options.

### EAM potential for solid copper

Call MD-Bench as follows:

```shell=
./MDBench-<TAG> -n 400 -f eam -e ./data/Cu_u3.eam
```

Two different EAM variants are available: `Cu_u3.eam` and `Cu_u6.eam`. The EAM
potential is currently only available for verletlist.

### Lennard-Jones potential for melted copper

The melted copper testcase has only 32000 atoms in the default configuration.
Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -n 400  -i ./data/copper_melting/input_lj_cu_one_atomtype_20x20x20.dmp
```

### EAM potential for melted copper

Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -n 400 -f eam -e ./data/Cu_u3.eam  -i ./data/copper_melting/input_eam_cu_one_atomtype_20x20x20.dmp
```

Two different EAM variants are available: `Cu_u3.eam` and `Cu_u6.eam`. The EAM
potential is currently only available for verletlist.

### Lennard-Jones potential for argon gas

Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -i ./data/argon/input.gro  -p ./data/argon/mdbench_params.conf
```

### GROMACS non-bonded kernel benchmark (clusterpair only)

`gmx nonbonded-benchmark` times the non-bonded kernels of GROMACS (NBNxM) on a box of
1000 water molecules with a fixed pair list. MD-Bench can run the same system with the
same physics (Lennard-Jones, Ewald or reaction-field electrostatics, intramolecular
exclusions, cut-off 1.0 nm) through its generic cluster-pair SIMD kernels, so forces and
kernel times can be compared. `tests/gromacs_nonbonded/` has the exact reference, the
GROMACS patch and the comparison script, see the README there.

Run the benchmark with MD-Bench (use one thread: half lists are not thread safe, and
GROMACS' benchmark also defaults to one):

```shell=bash
# 4xN kernel (SIMD=AVX2 and single precision give the 4x8 kernel of GROMACS' 4xM)
make OPT_SCHEME=clusterpair SIMD=AVX2 DATA_TYPE=SP CLUSTER_PAIR_KERNEL=4xN \
     COULOMB_RULE=ewald ACCURATE_RECIPROCAL=true
OMP_NUM_THREADS=1 ./MDBench-CP-4xN-GCC-X86-AVX2-SP-geometric-ewald-accrecip \
     -gmxbench -size 1 -half 1 --kernel-bench 500
```

`-size` replicates the box (a power of two, 3000 atoms times this value), `-half 1` uses
half neighbor lists, which is what GROMACS does, and `--kernel-bench` runs only the force
kernel that many times. The output ends with

```
KERNEL BENCHMARK
  Atoms:                 3000
  Cut-off:               1
  Neighbor lists:        half
  Cluster pairs in list: 45109 (4x8, 1443488 atom pairs)
  Iterations:            500
  Time per iteration:    1146.344 us
  Useful pairs per iter: 629741
  Computed pairs/us:     1259.210
  Useful pairs/us:       549.348
```

(The timings depend on the machine and vary by about 10% from run to run.)

The same with GROMACS (a `gmx` built with `tests/gromacs_nonbonded/gromacs-force-dump.patch`
adds the `-forces` option used to compare results):

```shell=bash
gmx nonbonded-benchmark -size 1 -nt 1 -simd 4xm -coulomb ewald -combrule geometric \
    -cutoff 1.0 -iter 500 -time -forces gmx_forces.txt
OMP_NUM_THREADS=1 ./MDBench-CP-4xN-GCC-X86-AVX2-SP-geometric-ewald-accrecip \
    -gmxbench -size 1 -half 1 --kernel-bench 500 --dump-forces mdbench_forces.txt
python3 tests/gromacs_nonbonded/compare_forces.py gmx_forces.txt mdbench_forces.txt
```

or build GROMACS, run both programs and compare forces and times in one go:

```shell=bash
GMX=/path/to/gmx PIN_CORE=3 tests/gromacs_nonbonded/run_comparison.sh 1 8
```

Other kernels: `CLUSTER_PAIR_KERNEL=2xNN` (compare with `-simd 2xmm`), `COULOMB_RULE=rf`
(compare with `-coulomb reaction-field`). `ACCURATE_RECIPROCAL=true` is needed for the
forces to agree with GROMACS (rms error about 1e-6 in single precision), without it
the kernels use the 12-bit hardware estimates and are a bit faster. `-gmxbench` is a
kernel benchmark: the molecules have no bonded terms, so integrating them for many steps
is not meaningful.

### Tabulated Lennard-Jones potential (verletlist only)

Uses cubic spline interpolation (GROMACS-style) instead of direct force computation.
Only available with `OPT_SCHEME=verletlist`. Select the indexing mode at build
time via `LJ_TABLE_INDEX=r` (default) or `LJ_TABLE_INDEX=rsq` (avoids sqrt,
load-bound):

```shell=bash
make OPT_SCHEME=verletlist LJ_TABLE_INDEX=rsq
./MDBench-<TAG> -f lj_table
```

## Citations

Rafael Ravedutti Lucio Machado, Jan Eitzinger, Jan Laukemann, Georg Hager,
Harald Köstler and Gerhard Wellein: MD-Bench: A performance-focused prototyping
harness for state-of-the-art short-range molecular dynamics algorithms. Future
Generation Computer Systems
([FGCS](https://www.sciencedirect.com/journal/future-generation-computer-systems)),
Volume 149, 2023, Pages 25-38, ISSN 0167-739X, DOI:
[https://doi.org/10.1016/j.future.2023.06.023](https://doi.org/10.1016/j.future.2023.06.023)

Rafael Ravedutti Lucio Machado, Jan Eitzinger, Harald Köstler, and Gerhard
Wellein: MD-Bench: A generic proxy-app toolbox for state-of-the-art molecular
dynamics algorithms. Accepted for [PPAM](https://ppam.edu.pl/) 2022, the 14th
International Conference on Parallel Processing and Applied Mathematics, Gdansk,
Poland, September 11-14, 2022. PPAM 2022 Best Paper Award. Preprint:
[arXiv:2207.13094](https://arxiv.org/abs/2207.13094), DOI:
[https://dl.acm.org/doi/10.1007/978-3-031-30442-2_24](https://dl.acm.org/doi/10.1007/978-3-031-30442-2_24)

## Credits

MD-Bench is developed by the Erlangen National High Performance Computing Center
([NHR@FAU](https://hpc.fau.de/)) at the University of Erlangen-Nürnberg.

## License

[LGPL-3.0](https://github.com/RRZE-HPC/MD-Bench/blob/master/LICENSE)
