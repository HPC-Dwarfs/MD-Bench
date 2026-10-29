# MD-Bench Agent Guide

This document provides essential information for agents working with the MD-Bench codebase, a molecular dynamics benchmarking tool.

## Project Overview

MD-Bench is a toolbox for performance engineering of short-range force calculation kernels used in molecular dynamics applications. It supports state-of-the-art algorithms from community codes like LAMMPS and GROMACS.

The codebase supports multiple optimization schemes:
- **verletlist**: Traditional Verlet neighbor list approach
- **clusterpair**: Modern cluster-pair approach based on GROMACS algorithms

The project supports both CPU and GPU implementations with various SIMD configurations.

## Code Structure

```
src/
├── clusterpair/     # Cluster-pair optimization scheme implementation
├── verletlist/     # Verlet list optimization scheme implementation
└── common/         # Shared utilities and common code
    └── simd/       # SIMD intrinsics for different architectures
```

Each optimization scheme has its own directory with similar file structures:
- `atom.*`: Atom data structures and management
- `force.*`: Force calculation kernels
- `neighbor.*`: Neighbor list construction
- `pbc.*`: Periodic boundary conditions
- `integrate.*`: Integration methods
- `main.*`: Main application entry points
- `main-stub.*`: Stub entry point for synthetic neighbor list benchmarking
- `stats.*`: Performance statistics
- `tracing.*`: Debug tracing functionality

## Build System

MD-Bench uses a sophisticated Make-based build system with automatic dependency generation.

### Essential Build Commands

```bash
# Build with default configuration
make

# Build with specific toolchain
make TOOLCHAIN=ICX

# Build with CUDA support (Nvidia)
make TOOLCHAIN=NVCC

# Build with HIP support (AMD); GPU_ARCH is required, no default
make TOOLCHAIN=HIPCC GPU_ARCH=gfx1030

# Clean build artifacts for current configuration
make clean

# Clean everything including binaries for all configurations
make cleanall

# Format source code
make format

# Generate assembly output
make asm

# Show compiler information
make info
```

### Binary Naming

Binary names follow `MDBench-<OPT_TAG>-<TOOL_TAG>-<DATA_TYPE>`, e.g.:
- `MDBench-VL-GCC-X86-AVX512-DP` (verlet-list, GCC, AVX512, double precision)
- `MDBench-CP-auto-GCC-X86-AVX512-DP` (cluster-pair, auto kernel selection, AVX512, DP)

Different configurations produce differently-tagged binaries that coexist in
the same directory.

### Configuration Options

Key configuration options in `config.mk`:

- `TOOLCHAIN`: Compiler toolchain (GCC/CLANG/ICC/ICX/NVCC/HIPCC)
- `OPT_SCHEME`: Algorithmic variant (verletlist/clusterpair)
- `SIMD`: SIMD instruction set (NONE/SSE/AVX/AVX2/AVX512/NEON/SVE/SVE2)
- `DATA_TYPE`: Precision (SP for single, DP for double)
- `ATOM_DATA_LAYOUT`: Data layout (AOS for array-of-structures, SOA for structure-of-arrays)
- `NBLIST_DATA_LAYOUT`: Neighbor-list data layout (auto/AOS/SOA); `auto` defaults to AOS for CPU and SOA for GPU
- `LJ_COMB_RULE`: LJ combination rule (`single` for global broadcast, `geometric` for per-type, `full` for full pair matrix)
- `COULOMB_RULE`: Electrostatics (`none`/`rf`/`ewald`), cluster-pair SIMD (4xN, 2xNN, 2xN) and GPU kernels; adds `-<rule>` to the binary tag
- `EXCLUSIONS`: Molecular exclusions (`auto`/`true`/`false`), `auto` = on when `COULOMB_RULE` is not `none`
- `ACCURATE_RECIPROCAL`: Newton-refine the 12-14 bit `simd_real_reciprocal`/rsqrt estimates used by the SIMD kernels, cluster-pair and verlet-list (`USE_SIMD_KERNEL`) (`src/common/simd_math.h`; default off, adds `-accrecip` to the tag)
- `ENABLE_MPI`: Enable MPI parallelization
- `ENABLE_OPENMP`: Enable OpenMP parallelization
- `GPU_ARCH`: Target GPU architecture for HIPCC (AMD) builds; **required, no default**. Common values: `gfx90a` (MI250X), `gfx940` (MI300A), `gfx942` (MI300X), `gfx1030` (RX 6800 XT)
- `LJ_TABLE_INDEX`: Indexing mode for tabulated LJ kernels (`r` or `rsq`); `rsq` avoids a sqrt and is load-bound
- `CLUSTER_PAIR_KERNEL`: Kernel variant for the clusterpair scheme (auto/4xN/2xNN/gpusimple)
- `ASM_SYNTAX`: Assembly syntax for generated `.s` files on X86 (`intel`/`att`), passed to the compiler as `-masm=`; ignored on ARM
- `USE_PINNED_MEMORY`: Use pinned (page-locked) host memory for buffers that are `cudaMemcpy`'d to/from the GPU every reneighbor step, for faster/more predictable transfers (default true; CUDA/HIP builds only)
- `SUPERCLUSTER_J_SLICES`: Number of thread-block slices the GPU grid is split into along each i-supercluster's j-neighbor list, for the super-clustering kernel (default 4)
- `ENABLE_LIKWID`: Enable LIKWID performance counter instrumentation

Example build commands:
```bash
# CPU build with AVX512 and double precision
make TOOLCHAIN=ICX SIMD=AVX512 DATA_TYPE=DP

# GPU build with CUDA
make TOOLCHAIN=NVCC

# Cluster-pair build with specific kernel
make TOOLCHAIN=ICX OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL=4xN
```

## Key Implementation Patterns

### Data Types

The codebase uses precision-aware data types:
- `MD_FLOAT`: Configurable precision floating point (float/double)
- `MD_FLOAT3`: 3-component vector with configurable precision
- `MD_FLOAT4`: 4-component vector with configurable precision

### SIMD Architecture

SIMD operations are implemented through architecture-specific headers in `src/common/simd/`. The implementation automatically selects the appropriate SIMD width based on the configuration:
- Single precision: 2x vector width of double precision
- Double precision: Base vector width (2 for NEON, 4 for AVX2, 8 for AVX512)

Several SVE2 double-precision half-neighbor horizontal operations (`simd_real_load_h_dual`, `simd_real_load_h_duplicate`, `simd_real_h_decr3`, etc.) are not yet implemented and will abort at runtime if called.

### Kernel Variants

For cluster-pair implementation, several kernel variants are available:
- Reference kernel (scalar, 4x4 clusters whatever the SIMD width, like the GROMACS
  plain-C kernel; also the kernel of `SIMD=NONE` builds)
- SIMD 4xN kernel (`force_lj_4xn.c`, SIMD width 2, 4 or 8 reals)
- SIMD 2xNN kernel (`force_lj_2xnn.c`, SIMD width 8 or 16 reals)
- SIMD 2xN kernel (`force_lj_2xn.c`, SIMD width 2 or 4 reals, for narrow registers such as NEON/SVE; runs on x86 only with AVX2 double precision)

`force.h` stops the build with an `#error` for any other width.
- CUDA/HIP GPU kernel (`gpusimple`)
- CUDA/HIP super-cluster kernel

Both schemes also have a tabulated LJ kernel (`force_lj_table.c`): cubic
spline interpolation (GROMACS-style) over a precomputed force/energy table.
Verlet-list additionally has a CUDA GPU variant (`force_lj_table_cuda.cu`);
cluster-pair's tabulated LJ is CPU-only.

The kernel selection is controlled through build-time options and affects the generated binary name.

Each scheme also builds a `-stub` binary variant (`main-stub.c`) for synthetic neighbor list benchmarking with configurable patterns (sequential, fixed, random).

### LJ Force Computation

Reference (scalar) force formula (`src/clusterpair/force_lj.c`,
`src/verletlist/force_lj.c`), applied when `rsq < cutforcesq`. The SIMD cluster-pair
kernels evaluate the same force as `(12*C12*r^-12 - 6*C6*r^-6) / r^2`, see below:
```c
MD_FLOAT sr2   = 1.0 / rsq;
MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
MD_FLOAT force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
```
Half neighbor-lists apply Newton's third law (each pair's force is added to
both atoms from a single computation).

Combination rules, selected at compile time via `-DLJ_COMB_RULE`
(`src/common/parameter.h`):
- `LJ_COMB_SINGLE` (0): single atom type, broadcasts global epsilon/sigma — fastest, no per-type lookup
- `LJ_COMB_GEOM` (1): per-type geometric combination — `epsilon_ij = sqrt(eps_i*eps_j)`, `sigma6_ij = sigma3_i*sigma3_j`
- `LJ_COMB_FULL` (2): full type-pair matrix lookup (gathers; the SIMD kernels read the `ljc6x6`/`ljc12x12` tables)

Cluster-pair kernel entry points (`force_lj.c` (reference), `force_lj_4xn.c`,
`force_lj_2xnn.c`, `force_lj_2xn.c`, `force_lj_cuda.cu` / `force_lj_cuda_sup.cu`):
`computeForceLJRef`,
`computeForceLJ4xnHalfNeigh`/`FullNeigh`, `computeForceLJ2xnnHalfNeigh`/`FullNeigh`,
`computeForceLJ2xnHalfNeigh`/`FullNeigh`, `computeForceLJCuda`, `computeForceLJCudaSup`.

### The SIMD cluster-pair kernels and their performance design

The CPU SIMD kernels live in `force_lj_4xn.c`, `force_lj_2xnn.c` and `force_lj_2xn.c`
(`computeForceLJ<kind>HalfNeigh`/`FullNeigh`); `force_lj.c` only has the scalar
reference kernel and the eight-shell ghost routine. The pair math is shared and written
once in `pair_rows.h` (`pairForces()`): a kernel keeps the i-cluster in "rows"
(`IRows`: broadcast atoms for 4xN/2xN, atom pairs in the register halves for 2xNN), loads
the j-cluster (`JData`) and only differs in how it loads, in the diagonal mask table and
in how it reduces forces. `kernel_math.h` pulls in the reciprocal/rsqrt helpers of `src/common/simd_math.h`
(`ACCURATE_RECIPROCAL` selects Newton-Raphson refinement). The math is GROMACS': one rsqrt
for 1/r and 1/r^2 (a plain reciprocal without electrostatics), the cut-off/diagonal mask
applied to 1/r only, LJ from per-atom `cl_sqrt6c6`/`cl_sqrt12c12` products or the
`ljc6x6`/`ljc12x12` type-pair tables, all phase by phase over the rows.
Things that matter for the speed and are easy to break:
- the GCC flags `-fschedule-insns -fsched-pressure -fno-unroll-loops` for the three kernel
  files (`KERNEL_CFLAGS` in `make/include_GCC.mk`, applied in the `Makefile`) are worth
  15%: GCC's default order of the dependency chains is poor. The functions in
  `pair_rows.h` are `always_inline` and must stay that way;
- padding atoms of a cluster are at `CLUSTER_PAD_COORD` (finite), `pairForces()` multiplies
  the zero force of a masked pair with its distance and would produce NaN with `INF`.
  Padding atoms of two clusters are at distance zero, so the kernels clamp the squared
  distance to `MIN_DISTANCE_SQ` after the cut-off check, like GROMACS
  (`c_nbnxmMinDistanceSquared`);
- j-cluster data of the entry 3 positions ahead in the visiting order (inner masked, then
  inner unmasked entries; `prefetchAhead()`) is prefetched (`KERNEL_PREFETCH=false`
  disables it, the binary gets the suffix `-nopf`). GROMACS does not prefetch, disable it
  for a like-for-like comparison;
- without MPI, `initForce()` selects the half-shell method (`param->method`) for half
  lists of every CPU kernel (before `initNeighbor()`, which reads it) and `setup()` sets
  `atom->mybox` to the simulation box: only the ghost clusters in the positive half space
  are kept, the kernels write the reaction forces into them (`cj < ... || param->method`)
  and `reverse()` folds them back. Every CPU kernel must scatter into ghost clusters.
- half lists write j-cluster forces that other threads own: the half-list kernels write
  through `threadForces()` (`thread_forces.h`), thread 0 into `cl_f` and the others into
  per-thread arrays that `threadForcesReduce()` sums into `cl_f` (GROMACS' per-thread
  output buffers). A single thread has no extra work.

### Electrostatics, exclusions and the GROMACS benchmark

`src/clusterpair/coulomb.h` holds the GROMACS-style SIMD Ewald correction
(analytical `pmeForceCorrection` rational approximation) and the constants of the reaction
field, used by `pairForces()` in all three SIMD kernels. The GPU kernels (`gpusimple` and
`supercluster`) evaluate each pair with `pairForceCuda()` (`src/clusterpair/pair_cuda.h`),
using the scalar copy of the correction in `src/common/coulomb_scalar.h`; they compare
`cl_mol` per pair, as their lists have no masked part. `cl_q`/`cl_mol` are only copied to the
device when `COULOMB_RULE`/`EXCLUSIONS` need them. `initForce()` aborts for the reference
kernel when `COULOMB_RULE`/`EXCLUSIONS` are enabled.
Charges are per atom type (`atom->q_per_type`, copied to `cl_q`). Exclusions are
molecule based: `atom->molsize` consecutive atom ids form a molecule, `cl_mol` holds the
molecule index per cluster slot and `cl_id` the original atom index (the atom array
itself is never permuted, only `updateSingleAtoms()` rewrites it from the clusters).
`buildNeighborCPU()` keeps cluster pairs that contain an excluded atom pair in the
masked part of the list (`clusterPairHasExclusion()`), the masked kernel loop compares
`cl_mol` per lane and only zeroes LJ and the direct `1/r^3` term: excluded pairs inside
the cut-off still get the exclusion correction, like in GROMACS' NBNxM kernels. The
unmasked loop has no exclusion work. Ghost clusters copy `cl_q`/`cl_mol`/`cl_id` in the
`ADDGHOST` macro (`pbc.c`). With MPI the atom ids are global (`initAtomIds()` runs before
the domain decomposition), the ghost exchange sends type and id per atom and
`unpackGhost()` derives the rest with `setClusterAtom()`, atom migration sends the id and
`--dump-forces` gathers on rank 0. This path is untested (no MPI on the dev machine).

`-gmxbench [-size N]` (`readAtomGmx()`, data in `src/common/gmxcoords.h`) builds the
water box of `gmx nonbonded-benchmark`. That benchmark passes C6/C12 where NBNxM expects
6*C6/12*C12, so its effective LJ potential uses C6/6 and C12/12; `readAtomGmx()` does the
same. `--kernel-bench N` times only the force kernel on the fixed lists and `--dump-forces
f` writes per-atom forces by original index. `tests/gromacs_nonbonded/` has the exact
float64 reference, the force comparison, the GROMACS patch and the comparison script.

## Core Data Structures

**Cluster-pair `Atom`** (`src/clusterpair/atom.h`) holds both flat per-atom
arrays (`x`/`y`/`z`, `vx`/`vy`/`vz`, `type`) and cluster-organized copies used
by the force kernels for cache locality: `cl_x`/`cl_v`/`cl_f`/`cl_t`, plus
per-atom geometric-combination-rule parameters in cluster layout
(`cl_sqrt_epsilon`, `cl_sigma3`). `iclusters`/`jclusters` are per-cluster
bounding-box + atom-count records (`Cluster`); `siclusters` are optional
super-cluster records (`SuperCluster`) used when `super_clustering` is
enabled. Coordinate offsets into `cl_x` are computed via macros
(`CL_X_INDEX`/`CL_Y_INDEX`/`CL_Z_INDEX`, `CI_BASE_INDEX`, etc. in
`src/clusterpair/force.h`) whose layout depends on the active kernel variant
(4xN vs 2xNN vs GPU) — always check which macros apply before changing
cluster indexing.

**Cluster-pair `Neighbor`** (`src/clusterpair/neighbor.h`) stores, per
i-cluster: `numneigh`/`numneigh_inner` (neighbor counts, inner = within the
pruned double-cutoff region), `neighbors` (j-cluster indices, AOS or SOA
depending on `NBLIST_DATA_LAYOUT`), and `neighbors_imask` (per-sub-cluster
interaction masks used by super-clustering and SIMD exclusion handling).

**`Parameter`** (`src/common/parameter.h`) carries all runtime configuration:
timestep (`dt`/`dtforce`), system size (`nx`/`ny`/`nz`), cutoffs
(`cutforce`/`skin`/`outer_skin`/`cutneigh`), LJ parameters
(`epsilon`/`sigma`/`sigma6`, or `epsilon_per_type`/`sigma_per_type` plus
`ntypes` for multi-type systems), `LJ_COMB_RULE`, `half_neigh`,
`reneigh_every`, and GPU-specific flags like `super_clustering`. Loaded from
a `key value` text file (`readParameter`/`readTypesFile` in
`src/common/parameter.c`); the types file is two columns (epsilon, sigma).

## Key Files for Development

| File | Purpose |
|------|---------|
| `src/common/parameter.h`, `parameter.c` | Runtime parameter struct and file parsing |
| `src/clusterpair/atom.h` | Cluster-pair atom/cluster data structures |
| `src/verletlist/atom.h` | Verlet-list atom structure (simpler) |
| `src/clusterpair/neighbor.h` | Cluster-pair neighbor list structure |
| `src/clusterpair/force.h` | Kernel macro definitions and cluster indexing |
| `src/clusterpair/force_lj.c` | Cluster-pair scalar reference kernel, eight-shell ghost routine |
| `src/clusterpair/force_lj_4xn.c`, `force_lj_2xnn.c`, `force_lj_2xn.c` | SIMD kernels (half/full list) |
| `src/clusterpair/pair_rows.h` | pair math of the SIMD kernels: LJ, Coulomb, exclusions |
| `src/common/simd_math.h`, `src/clusterpair/kernel_math.h`, `coulomb.h` | accurate reciprocal/rsqrt; Ewald (`pmeForceCorrectionSimd`) |
| `src/clusterpair/thread_forces.h` | per-thread force arrays of the half-list kernels |
| `tests/gromacs_nonbonded/` | GROMACS benchmark comparison: reference, patch, scripts |
| `src/clusterpair/force_lj_cuda.cu`, `force_lj_cuda_sup.cu` | Cluster-pair CUDA/HIP LJ kernels (plain and super-clustering) |
| `src/verletlist/force_lj.c` | Verlet-list LJ force kernels |
| `src/common/ljtable.c` | Tabulated LJ force lookup |
| `src/common/device.c`, `device.h`, `deviceGPU.cu` | CUDA/HIP allocation, memcpy, and pinned-memory helpers |
| `src/clusterpair/main.c`, `src/verletlist/main.c` | Main program entry (argument parsing, simulation loop) |
| `config.mk` | Build configuration variables |
| `Makefile` | Build rules (pattern rules, dependency generation, test targets) |
| `.clang-format` | Code formatting rules (90-char limit, WebKit style) |

## Testing and Execution

### Running the Application

Execute with default parameters:
```bash
./MDBench-<TAG>
```

Common command line arguments:
- `-n <steps>`: Number of timesteps (default 200)
- `-nx/-ny/-nz <int>`: System dimensions (default 32x32x32)
- `-f <field>`: Force field (`lj` for Lennard-Jones, `eam` for EAM, `lj_table` for tabulated LJ)
- `-half <0|1>`: Use half (1) or full (0) neighbor lists
- `-r <radius>`: Cutoff radius (default 2.5)
- `-s <skin>`: Skin for Verlet buffer (default 0.3)
- `-p <file>`: Parameter file; `outer_skin` can be set here to enable double-cutoff pruning
- `-gmxbench [-size N]`, `--kernel-bench N`, `--dump-forces <file>`: the water system of the GROMACS non-bonded benchmark, a kernel-only timing loop and a per-atom force dump (cluster-pair only, see below)

### Test Cases

Several built-in test cases are available:
1. Lennard-Jones potential for solid copper (default)
2. EAM potential for solid copper: `-f eam -e ./data/Cu_u3.eam`
3. Lennard-Jones for melted copper with input file
4. Argon gas simulation with parameter files
5. Water box of the GROMACS non-bonded benchmark: `./MDBench-<TAG> -gmxbench -size 1 -half 1 --kernel-bench 500`
   (build with `COULOMB_RULE=ewald ACCURATE_RECIPROCAL=true`; `OMP_NUM_THREADS=1` like GROMACS' default)

### Performance Measurement

The application outputs performance metrics including:
- Time per step
- Performance in million atom-updates per second (MA/s)
- Memory bandwidth usage
- Hardware performance counters (when enabled)

### Automated Tests

Unit tests compile into one binary per configuration, `<benchmark binary>-tests` (not
part of the default `make` target), covering parameter parsing, atom/data layout, force
kernels, neighbor list construction, integration, box geometry, and thermodynamics. In
single precision `ASSERT_NEAR` uses a relative tolerance of at least 1e-5:
```bash
make OPT_SCHEME=clusterpair mdbench-tests
./MDBench-CP-auto-GCC-X86-AVX512-SP-geometric-tests
```

Regression tests are shell scripts under `tests/` (e.g.
`sim_argon_regression.sh`, `regression_scheme_equiv.sh`,
`test_simd_vs_scalar.sh`, `test_data_layout.sh`) that build one or more
binary configurations and compare their output. Two of them compare per-atom forces:
`test_halfshell.sh` (half vs. full lists for 4xN, 2xNN and 2xN, argon and water) and
`test_coulomb_exclusions.sh` (reaction-field, Ewald and exclusions of the three SIMD
kernels against the exact float64 reference in `tests/gromacs_nonbonded/`; both need
python3 with numpy, the second also scipy, and skip otherwise). `test_gpu_coulomb.sh` is
the same check for the GPU kernels; it needs nvcc and a GPU and is not part of `make
test`. With the NVIDIA HPC SDK module, nvcc and the CUDA runtime are not on the default
paths: add `$NVHPC_ROOT/Linux_x86_64/<version>/cuda/<cuda version>/bin` to `PATH` and
its `lib64` to `LD_LIBRARY_PATH`. `SUPERCLUSTER_INVERSE_THREAD_MAPPING` is not part of the
binary tag, so switching it needs a clean build. Run everything with:
```bash
make test                                  # unit tests + all regression scripts
make SIMD=NONE test                        # test a specific configuration
bash tests/sim_argon_regression.sh ./MDBench-VL-GCC-X86-AVX512-DP
```
Apart from `test_gpu_coulomb.sh`, which has to be run by hand on a GPU node, there is no
automated coverage for the cluster-pair GPU/super-clustering path
(`CLUSTER_PAIR_KERNEL=gpusimple`/`supercluster`).

Comparing with a real GROMACS needs a build of `gmx` with
`tests/gromacs_nonbonded/gromacs-force-dump.patch` and is done by
`tests/gromacs_nonbonded/run_comparison.sh` (not part of `make test`). Performance work
on the SIMD kernels should be measured with `perf stat` (cycles, instructions, FP pipe
and dispatch-stall events) on `--kernel-bench` runs, not with wall-clock times, which
vary with the clock speed; `tests/gromacs_nonbonded/README.md` lists what was tried.

Simulation test data lives in `data/`: `data/argon/input.gro` +
`mdbench_params*.conf` (argon gas, GROMACS `.gro` format) and
`data/copper_melting/` (copper FCC system).

## Important Gotchas

1. **Build Tags**: Each configuration produces a uniquely named binary with a tag reflecting the build options. Multiple configurations can coexist.

2. **CUDA Limitations**: CUDA builds only support GPU kernels; CPU kernels are disabled when TOOLCHAIN=NVCC.

3. **HIP Support**: Both verlet-list and cluster-pair support HIP (HIPCC) builds. HIP sources are auto-generated via `hipify-perl` from the same `.cu` sources used for CUDA at build time; those sources use `#ifdef __HIP_PLATFORM_AMD__` where AMD-specific behavior is needed. `GPU_ARCH` is always required and has no default.

4. **GPU PBC**: `updatePbcCUDA` exists only for verlet-list (`src/verletlist/pbcCuda.cu`). In cluster-pair CUDA builds, PBC still runs on the CPU.

5. **Double Cutoff Pruning**: Enabled by setting `outer_skin > 0` in a parameter file. The outer cutoff = `cutforce + skin + outer_skin`. GPU pruning kernels run device-resident for both verlet-list and cluster-pair CUDA builds. The inner-section tail padding for `CLUSTER_N < VECTOR_WIDTH` (2xNN) is correct but may visit extra j-clusters (filtered by per-pair cutforcesq).

6. **SIMD Dependencies**: SIMD kernels require specific hardware support. The reference kernel is used automatically when SIMD is not available.

7. **Data Layout Impact**: AOS vs SOA layout affects both performance and memory access patterns significantly.

8. **Neighbor List Construction**: The choice between full and half neighbor lists affects both correctness and performance.

9. **Precision Effects**: Single vs double precision not only affects accuracy but also vectorization efficiency.

10. **Header Dependencies**: The force.h header contains complex macro definitions that control kernel generation.

11. **Feature Parity Gaps**: Tabulated LJ CPU kernels exist for both schemes (`force_lj_table.c`), but the CUDA-accelerated variant (`force_lj_table_cuda.cu`) exists only for verlet-list. GPU PBC exists only for verlet-list (`pbcCuda.cu`); cluster-pair CUDA/HIP builds still run PBC on the CPU. EAM exists for both schemes but has no GPU kernel in either — it always runs on CPU regardless of `TOOLCHAIN`. Electrostatics/exclusions exist only for the SIMD CPU and GPU cluster-pair kernels, not for the reference kernel and verlet-list.

12. **Multiple Kernels Per Change**: Don't assume a single force kernel; a correctness or feature change in cluster-pair may need to be applied across reference, 4xN, 2xNN, 2xN, CUDA, and super-clustering GPU variants to stay consistent.

13. **Per-Type Parameters**: `LJ_COMB_GEOM` requires per-type `sqrt(epsilon)` and `sigma^3` (both flat, per-type, and cluster-layout copies); `LJ_COMB_FULL` additionally needs full type-pair epsilon/sigma6 matrices, which are not SIMD-optimized.

## Development Process

1. **Modify source files**: Changes are automatically detected through dependency tracking
2. **Build**: Use appropriate make command with configuration options
3. **Test**: Run with appropriate test case parameters
4. **Format**: Use `make format` to maintain code style consistency
5. **Clean**: Use `make clean` or `make cleanall` as needed

## Code Style and Conventions

- C11 standard with compiler-specific extensions for SIMD and GPU
- Function names use camelCase
- Constants and macros use UPPER_CASE
- Struct typedefs use PascalCase (`Atom`, `Parameter`, `Neighbor`, ...) — no `_t` suffix
- File-local functions use `static` keyword
- Extensive use of preprocessor conditionals for feature selection
- Heavy use of macros for code generation in performance-critical kernels

The project uses `.clang-format` for code formatting consistency.
