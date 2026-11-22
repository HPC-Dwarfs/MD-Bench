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
`config.mk`, for starters on a X86 based system the defaults are fine. `config.mk` only
holds the options; `make/options.mk` turns them into compiler defines and the binary
name.

To build the binary call, don't forget to load the compiler module on the
NHR@FAU clusters (e.g. `module load intel`):

```shell=
make
```

While the Makefile works with any version of GNU make, some features require GNU
make > v4.

You can run MD-Bench without any arguments:

```shell=
./MDBench-VL-ICX-AVX2-DP-omp
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

Multiple configurations can be build at the same time. The binary is named
`./MDBench-<build tag>` and the intermediate build results are located in a
`./build/build-<build tag>/` directory. The tag holds the options that make a binary a
different benchmark:

```
MDBench-<scheme>-<toolchain>[-<SIMD>][-CMP]-<SP|DP>[-<Coulomb>][-omp][-mpi]
```

`<scheme>` is `VL` (with `-SIMD` and `-NBSIMD` for `USE_SIMD_KERNEL` and
`USE_SIMD_NEIGHBOR`) or `CP-<CLUSTER_PAIR_KERNEL>`, where the C reference kernel of the CPU
is `CP-ref`. Other options (data layouts, tracers, statistics, GPU settings, ...) are not
part of the name, builds that differ only in them write the same binary. Every build
directory records the compiler and all flags in `config.stamp`, a build with other options
recompiles everything, so a binary always matches the options of its last build. Use a
different `PREFIX_DIR` to keep such builds side by side. `make -s <options> print-target`
prints the path of the binary of a configuration, for scripts.

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
- `make print-target`: Print the path of the binary of the configuration.
- `make asm`: Generate assembly output of all source files. The assembly files
are placed in the intermediate build directory.
- `make format`: Reformat all source files with `clang-format` using the format
specification in `.clang-format`.

### Build time options

- `TOOLCHAIN`: Determines which toolchain makefile is included
(`GCC`/`CLANG`/`ICC`/`ICX`/`ONEAPI`/`NVCC`/`HIPCC`)
- `SIMD`: Instruction set of the build, it also selects the architecture (there is no
`ISA` option anymore). X86: `X86_SCALAR`/`SSE`/`AVX`/`AVX_FMA`/`AVX2`/`AVX512`; ARM:
`ARM_SCALAR`/`NEON`/`SVE`/`SVE2`. `X86_SCALAR`/`ARM_SCALAR` is plain scalar code without
SIMD extensions (the C reference kernels), `NONE` the scalar code of the build machine.
`AVX` builds without FMA, `AVX_FMA` with it (both without AVX2 instructions). The
cluster-pair kernels and the verletlist SIMD kernels have no SSE version, make stops for
them with `SIMD=SSE`. GPU builds (`NVCC`/`HIPCC`) run no CPU SIMD code and have no SIMD
part in the binary name.
- `OPT_SCHEME`: Algorithmic variant (`verletlist` or `clusterpair`), different
source directories and main routines are used
- `ENABLE_LIKWID`: Turn on LIKWID instrumentation, the LIKWID library has to be available
- `ENABLE_OPENMP`: Enable OpenMP thread-level parallelization (true or false, adds
`-omp` to the binary name)
- `ENABLE_MPI`: Turn on the MPI parallel version of the code (adds `-mpi` to the binary
name)
- `DATA_TYPE`: Switch between single precision (`SP`) and double precision (`DP`) floating
point. This is controlled by defines.
- `ATOM_DATA_LAYOUT`: Switch between array-of-structure (`AOS`) and structure-of-array
(`SOA`) layout for atom positions and forces. Tradeoff between better cache
utilisation and easier SIMD vectorization.
- `NBLIST_DATA_LAYOUT`: Neighbor-list data layout (`auto`/`AOS`/`SOA`). `auto` defaults
to `AOS` for CPU builds and `SOA` for GPU builds.
- `COULOMB_RULE`: Electrostatics (`none` (default)/`rf`/`ewald`/`ewald_table`) for all
cluster-pair kernels (SIMD, reference and GPU) and for the verletlist plain C and GPU
kernels. `rf` is a reaction-field with infinite dielectric, `ewald` the real-space Ewald
part with the analytical force correction used by GROMACS, `ewald_table` the same with
GROMACS' tabulated correction (`gmx nonbonded-benchmark -table`, CPU kernels only). The
verletlist explicit SIMD kernels (`USE_SIMD_KERNEL=true`), EAM and the tabulated LJ
kernels do not support it. Atom charges are given per atom type.
- `DEBUG`: Enable additional debug output
- `USE_REFERENCE_KERNEL`: Enforce usage of the C reference implementation of the
cluster-pair scheme for validation (automatically set for GPU and scalar builds,
binary name `CP-ref` on the CPU). It uses 4x4 clusters whatever the SIMD width.
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
- `SUPERCLUSTER_J_SLICES` (default: 0, auto): Number of thread-block slices the GPU
grid is split into along the j-neighbor lists of the optimized GPU kernels
(gpusimple, supercluster and verletlist). 1 means each i-(super)cluster's neighbors
are handled by a single block; >1 splits the neighbor list across multiple concurrent
blocks (more parallelism/occupancy) that combine partial forces with
atomicAdd (more atomic contention); 0 chooses the slices per force call from the
occupancy and the list sizes (default).
- `XTC_OUTPUT`: Enable XTC trajectory output (GROMACS file format)
- `MEM_TRACER`: Trace memory addresses for cache simulation (true or false)
- `INDEX_TRACER`: Trace indexes and distances for gather-MD analysis (true or false)
- `COMPUTE_STATS`: Compute and print additional runtime statistics. The C kernels (the
verletlist kernel without `USE_SIMD_KERNEL` and the cluster-pair reference kernel) also
count the pairs inside and outside of the cut-off.

The former build options `LJ_COMB_RULE`, `ACCURATE_RECIPROCAL`, `KERNEL_PREFETCH` and
`SORT_ATOMS` are the runtime options `--lj-comb`, `--accurate-recip`, `--no-prefetch` and
`--sort-every` now, `ISA` follows from `SIMD` and `EXCLUSIONS` from the system (see
below); make stops with a hint when they are given on the command line or in `config.mk`.

Molecular exclusions are on when the system has molecules (`atom->molsize > 1`, the
water of `-gmxbench`): atoms of the same molecule (consecutive atom ids) do not interact
via LJ or the direct Coulomb term, but get the exclusion correction of the selected
`COULOMB_RULE`. The kernels have a compiled variant for systems with molecules, so systems
without them pay nothing.

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
32x32x32 unit cells. Both schemes accept the same options (`-h` lists them), options
that only one scheme implements stop the other one with an error, and so does an unknown
option.

Input and system:

- `-p / --param <file>`: file to read parameters from (can be specified more
than once, `--params` is accepted as well). Default initialization sets parameters for
default LJ testcase. Supports `outer_skin = <value>` to enable double-cutoff pruning
(outer cutoff = `cutforce + skin + outer_skin`; default 0.0 disables it).
- `-f <string>`: force field (`lj`, `eam`, `lj_table`), default `lj`. For
anything different than `lj` you also need to provide a specific parameter file.
`lj_table` uses tabulated cubic Hermite spline interpolation. It is available for
`OPT_SCHEME=verletlist`, and for `OPT_SCHEME=clusterpair` only when built with the
reference kernel (`USE_REFERENCE_KERNEL=true`).
- `-i <file>`:  input file with atom positions. MD-Bench supports
Brookhaven protein data bank (.pdb), GROMACS GROMOS87 (.gro), and LAMMPS dump
(.dmp) file formats
- `-e <file>`:  input file for EAM parameters
- `-nx/-ny/-nz <int>`:  set linear dimension of systembox in x/y/z direction
(default 32 in every dimension)
- `-m / --mass <real>`: atom mass
- `--ntypes <int>`: number of atom types
- `-setup <int>`: 1 (default) creates the FCC lattice of the argon/copper setup, 0 skips
it
- `-gmxbench`:  use the water system of the GROMACS non-bonded benchmark
(`src/gromacs/nbnxm/benchmark`) instead of the default testcase: 3000 SPC/E-like
atoms with charges and intramolecular exclusions, cut-off 1.0 and no skin unless
`-r`/`-s` are given. `-size <int>` replicates the box (power of two) as `gmx
nonbonded-benchmark -size` does. Available in both schemes. Build with `COULOMB_RULE=ewald`
(or `rf`) and `CLUSTER_PAIR_KERNEL=4xN` to compare with GROMACS. It is a force kernel benchmark like GROMACS', so it needs `--kernel-bench` or `-n
0`.
- `-w <file>`:  write input atoms to file (verletlist)

Simulation and neighbor lists:

- `-n / --nsteps <int>`:  set number of timesteps for simulation (default 200)
- `-r / --radius <real>`:   set cutoff radius (default 2.5)
- `-s / --skin <real>`:   set skin (verlet buffer, default 0.3)
- `--outer-skin <real>`: additive outer skin, > 0 enables the dual pair list
- `--prune-every <int>`: steps between prunings of the inner list (default 5)
- `--reneigh-every <int>`: steps between neighbor-list builds (default 20)
- `--displacement-reneigh`: rebuild the lists when an atom moved more than skin/2 (use
`--reneigh-every 0` to disable the fixed interval)
- `--sort-every <int>`: sort the atoms spatially every `<int>` steps, followed by a
neighbor-list build; 0 disables it (verletlist; the cluster-pair scheme sorts
the atoms into clusters at every build). Default: off on the CPU, while the GPU
verletlist without MPI sorts the atoms on the device at every neighbor-list build.
Parameter file key `resort_every`.
- `-half <int>`: use half (1) or full (0) neighbor lists (default 0 - full
neighbor list). The cluster-pair SIMD kernels without MPI keep only half of
the periodic ghost clusters with half lists, so a half list has about half the entries of
a full list. With several OpenMP threads the half-list kernels accumulate into per-thread
force arrays that are summed afterwards
- `-method <int>`: ghost interactions with MPI: 0 full shell, 1 half shell, 2 eighth
shell, 3 half stencil (verletlist)
- `-bal <int>`: load balancing with MPI: 0 none, 1 RCB, 2 mean time RCB, 3 staggered

Force kernel:

- `--lj-comb <string>`: Lennard-Jones combination rule. `full` (default) reads the
type-pair tables (gathers), `geometric` combines per-type parameters (sqrt(eps_i eps_j),
sigma3_i sigma3_j), `single` uses one set of global parameters (fastest, no type lookup;
all atom types must have the global epsilon and sigma), `auto` picks `single` when all
types are alike and `geometric` otherwise. The type-pair tables are the geometric
combination of the per-type parameters, so the rules differ in speed, not in the forces.
The neighbor-list builds check the per-type cut-off for `geometric`/`full` and the global
one for `single`.
Parameter file key `lj_comb_rule`.
- `-energy`:  also compute the LJ and Coulomb energies in the force kernels (GROMACS'
energy kernels, `gmx nonbonded-benchmark -energy`), reported by `--kernel-bench` and
`--dump-forces`. All LJ kernels except the verletlist explicit SIMD kernels.
- `-halflj`:  LJ only for the first half of the i-clusters that have at most half of
their atoms with LJ, whose LJ atoms are sorted to the front (`gmx nonbonded-benchmark
-halflj`). CPU cluster-pair kernels only.
- `--no-accurate-recip`: use the hardware reciprocal and reciprocal square root estimates
(12-14 bits) of the SIMD kernels (cluster-pair and verletlist with `USE_SIMD_KERNEL`)
without refining them. By default they are refined with Newton-Raphson steps, as GROMACS
does, so the forces agree with a code that evaluates 1/r^2 to full precision and between
the SIMD widths (`--accurate-recip` selects the default explicitly). On the GPU it
selects the fast division and rsqrt (as `-use_fast_math`) instead of the IEEE ones, per
kernel, so both versions are in the same binary. Parameter file key
`accurate_reciprocal` (0 or 1).
- `-ref`: GPU builds: use the reference force kernels (one thread per atom pair, no
packing, j-slicing or warp reductions, all forces with atomicAdd) instead of the
optimized ones. Parameter file key `ref_kernel`.
- `--no-prefetch`: no software prefetch of the j-cluster data in the CPU cluster-pair
SIMD kernels. GROMACS does not prefetch, use it for a like-for-like comparison.
Parameter file key `kernel_prefetch`.
- `--kernel-bench <int>`:  after setup only run the force kernel this many times on the
fixed neighbor lists and report the time per call and the useful pair rate, as the
GROMACS benchmark does
- `--dump-forces <file>`:  write the forces of a single kernel pass (`atom id fx fy
fz`, original atom order) to a text file, with `-energy` also the energies

The force kernels and the neighbor-list builds of both schemes are templates compiled
once per variant, i.e. per combination of `--lj-comb`, `--accurate-recip`, `--no-prefetch`
and molecular exclusions as far as a kernel has them (`src/common/kernel_variant.h`; the
variants are the `<file>.<variant>.o` objects in the build directory, the GPU kernels are
C++ templates), so these options select a variant and add no work to the kernel.

Output:

- `--freq <real>`:  processor frequency (GHz), used to calculate cycle metrics
(default 2.4)
- `--vtk <file>`:    VTK output file for visualization
- `--xtc <file>`:    XTC trajectory (clusterpair, build with `XTC_OUTPUT=true`)
- `-v / --verbose`: per-routine timing breakdown in the performance report

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
GROMACS patch and the comparison script.

Run the benchmark with MD-Bench (one thread, the default of GROMACS' benchmark):

```shell=bash
# 4xN kernel (SIMD=AVX2 and single precision give the 4x8 kernel of GROMACS' 4xM)
make OPT_SCHEME=clusterpair SIMD=AVX2 DATA_TYPE=SP CLUSTER_PAIR_KERNEL=4xN \
     COULOMB_RULE=ewald
OMP_NUM_THREADS=1 ./MDBench-CP-4xN-GCC-AVX2-SP-ewald-omp \
     -gmxbench -size 1 -half 1 --kernel-bench 500 --accurate-recip --lj-comb geometric
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
OMP_NUM_THREADS=1 ./MDBench-CP-4xN-GCC-AVX2-SP-ewald-omp \
    -gmxbench -size 1 -half 1 --kernel-bench 500 --accurate-recip --lj-comb geometric \
    --dump-forces mdbench_forces.txt
python3 tests/gromacs_nonbonded/compare_forces.py gmx_forces.txt mdbench_forces.txt
```

or build GROMACS, run both programs and compare forces and times in one go:

```shell=bash
GMX=/path/to/gmx PIN_CORE=3 tests/gromacs_nonbonded/run_comparison.sh 1 8
```

Other kernels: `CLUSTER_PAIR_KERNEL=2xNN` (compare with `-simd 2xmm`; `run_comparison.sh`
does with `SIMD=AVX512` and a `gmx` built with `GMX_SIMD=AVX_512`), `COULOMB_RULE=rf`
(compare with `-coulomb reaction-field`), `COULOMB_RULE=ewald_table` (compare with
`-table`). The options `-energy` and `-halflj` are the same in both programs.
The refined reciprocals (default) are needed for the forces to agree with GROMACS (rms
error about 1e-6 in single precision, energies about 1e-6 relative), with
`--no-accurate-recip` the kernels use the 12-bit hardware estimates and are a bit faster.
`--no-prefetch` turns off the j-cluster prefetch, which GROMACS does not have. `--lj-comb geometric` is GROMACS' `-combrule geometric`
(the default `full` reads type-pair tables and is slower). `-gmxbench` is a kernel benchmark: the
molecules have no bonded terms or constraints, so MD-Bench refuses to integrate them
(use `--kernel-bench` or `-n 0`).

The energies follow the GROMACS benchmark settings: no potential shift for LJ and Ewald,
the `c_rf` shift for reaction-field, the exclusion correction for excluded pairs inside
the cut-off and the self energy of every charge.

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
