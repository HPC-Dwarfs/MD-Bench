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
- Reference kernel (scalar)
- SIMD 4xN kernel
- SIMD 2xNN kernel
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
`src/verletlist/force_lj.c`), applied when `rsq < cutforcesq`:
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
- `LJ_COMB_FULL` (2): full type-pair matrix lookup — not SIMD-optimized

Cluster-pair kernel entry points (`src/clusterpair/force_lj.c` /
`force_lj_cuda.cu` / `force_lj_cuda_sup.cu`): `computeForceLJRef`,
`computeForceLJ4xnHalfNeigh`/`FullNeigh`, `computeForceLJ2xnnHalfNeigh`/`FullNeigh`,
`computeForceLJ2xnHalfNeigh`/`FullNeigh`, `computeForceLJCuda`, `computeForceLJCudaSup`.

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
| `src/clusterpair/force_lj.c` | Cluster-pair LJ force kernels (reference, SIMD) |
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

### Test Cases

Several built-in test cases are available:
1. Lennard-Jones potential for solid copper (default)
2. EAM potential for solid copper: `-f eam -e ./data/Cu_u3.eam`
3. Lennard-Jones for melted copper with input file
4. Argon gas simulation with parameter files

### Performance Measurement

The application outputs performance metrics including:
- Time per step
- Performance in million atom-updates per second (MA/s)
- Memory bandwidth usage
- Hardware performance counters (when enabled)

### Automated Tests

Unit tests compile into a single `mdbench-tests` binary (not part of the
default `make` target) covering parameter parsing, atom/data layout, force
kernels, neighbor list construction, integration, box geometry, and
thermodynamics:
```bash
make mdbench-tests
./mdbench-tests
```

Regression tests are shell scripts under `tests/` (e.g.
`sim_argon_regression.sh`, `regression_scheme_equiv.sh`,
`test_simd_vs_scalar.sh`, `test_data_layout.sh`) that build one or more
binary configurations and compare their output. Run everything with:
```bash
make test                                  # unit tests + all regression scripts
make SIMD=NONE test                        # test a specific configuration
bash tests/sim_argon_regression.sh ./MDBench-VL-GCC-X86-AVX512-DP
```
There is no automated coverage for the cluster-pair GPU/super-clustering path
(`CLUSTER_PAIR_KERNEL=gpusimple`/`supercluster`) — those require manual
verification on real GPU hardware.

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

11. **Feature Parity Gaps**: Tabulated LJ CPU kernels exist for both schemes (`force_lj_table.c`), but the CUDA-accelerated variant (`force_lj_table_cuda.cu`) exists only for verlet-list. GPU PBC exists only for verlet-list (`pbcCuda.cu`); cluster-pair CUDA/HIP builds still run PBC on the CPU. EAM exists for both schemes but has no GPU kernel in either — it always runs on CPU regardless of `TOOLCHAIN`.

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
