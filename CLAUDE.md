# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

MD-Bench is a molecular dynamics simulation kernel benchmark suite with highly configurable build system.

### Build Configuration

All configuration is in `config.mk`. Key variables:

**Optimization scheme:**
- `OPT_SCHEME` (default: verletlist): `verletlist` or `clusterpair`. Changes the algorithm and source directory (`src/verletlist/` vs `src/clusterpair/`).

**Hardware/compiler:**
- `TOOLCHAIN` (default: GCC): GCC/CLANG/ICC/ICX/ONEAPI/NVCC/HIPCC
- `ISA` (default: X86): X86 or ARM
- `SIMD` (default: AVX512): X86: NONE/SSE/AVX/AVX_FMA/AVX2/AVX512; ARM: NONE/NEON/SVE/SVE2

**Precision & layout:**
- `DATA_TYPE` (default: SP): SP (single) or DP (double precision)
- `ATOM_DATA_LAYOUT` (default: AOS): AOS (array-of-structures, interleaved x/y/z) or SOA (struct-of-arrays)
- `NBLIST_DATA_LAYOUT` (default: auto): AOS (atom-major) or SOA (neighbor-major); auto selects based on target (AOS for CPU, SOA for GPU)

**Kernel variants (cluster-pair only):**
- `CLUSTER_PAIR_KERNEL` (default: auto): auto/4xN/2xNN/gpusimple
- `SUPERCLUSTER_DATA_LAYOUT` (default: AOS3): AOS3/AOS4/SOA (for super-clustering CUDA kernels)

**Other features:**
- `ENABLE_LIKWID` (default: false): Enable LIKWID performance counter instrumentation
- `ENABLE_OPENMP` (default: true): Enable OpenMP parallelization
- `ENABLE_MPI` (default: false): Enable MPI distributed parallelization
- `LJ_COMB_RULE` (default: geometric): single/geometric/full (LJ parameter combination rule)
- `LJ_TABLE_INDEX` (default: r): r or rsq (tabulated force lookup grid spacing)
- `USE_REFERENCE_KERNEL` (default: false): Use reference scalar kernel instead of SIMD
- `USE_SIMD_KERNEL` (default: false): Force SIMD intrinsic kernels (when available)
- `USE_SIMD_NEIGHBOR` (default: false): Use SIMD intrinsics to build the Verlet-list neighbor lists. Independent of `USE_SIMD_KERNEL` (which only affects the force kernel). Supported for AVX2/AVX512/NEON/SVE/SVE2 (not SSE/AVX); not implemented for `NBLIST_DATA_LAYOUT=CSR`. Note: the NEON and SVE/SVE2 code paths have not been compiled or run on real ARM hardware/toolchain; only AVX2/AVX512 have been build- and regression-verified.
- `SORT_ATOMS` (default: false): Sort atoms by frequency
- `MEM_TRACER`/`INDEX_TRACER`/`COMPUTE_STATS`: Enable performance tracing/statistics

### Binary Naming

Binary names are: `MDBench-<OPT_TAG>-<TOOL_TAG>-<DATA_TYPE>`

Examples:
- `MDBench-VL-GCC-X86-AVX512-DP` (verletlist with GCC, AVX512, double precision)
- `MDBench-CP-auto-GCC-X86-AVX512-DP` (clusterpair with auto kernel selection, AVX512, DP)

### Build Commands

```bash
# Build with current config.mk defaults
make

# Build with specific configuration (override config.mk)
make TOOLCHAIN=GCC OPT_SCHEME=clusterpair SIMD=AVX512 DATA_TYPE=DP

# Clean current build
make clean

# Clean all builds and generated files
make cleanall

# Generate assembly for current config
make asm

# Reformat all source files with clang-format (uses .clang-format spec)
make format

# Output compiler version (useful for logging)
make info

# Run full test suite (unit tests + regression tests)
make test
```

### Build Directory Structure

- `./src/common/`: Shared code (parameter parsing, atom data, box, thermodynamics, I/O, LJ table)
- `./src/verletlist/`: Verlet list scheme implementation
- `./src/clusterpair/`: Cluster-pair scheme (GROMACS nbnxn-style, CPU + GPU variants)
- `./make/`: Toolchain-specific build rules (include_GCC.mk, include_NVCC.mk, etc.)
- `./build/build-<TAG>/`: Intermediate object files and dependency files for each configuration
- `.clang-format`: Clang formatting specification (WebKit-based style, 90-char column limit)

## Architecture Overview

### Two Optimization Schemes

**Verlet List (`src/verletlist/`)**
- Traditional atom-centric neighbor list: for each atom, store list of neighbors within cutoff
- Simple data layout: per-atom position (x, y, z), velocity, force, type arrays
- Scalar or SIMD kernels available
- Simpler implementation, easier for small-scale parallelism

**Cluster Pair (`src/clusterpair/`)**
- GROMACS nbnxn-style algorithm: organizes atoms into small clusters, builds neighbor lists of clusters
- Two cluster dimensions: M (typically 4) and N (VECTOR_WIDTH, e.g., 8 for AVX512)
- Supports two kernel layouts: 4xN (calculates more pairs, better scheduling) and 2xNN (fewer pairs, less work)
- Data stored in cluster format (cl_x, cl_v, cl_f) for cache locality
- GPU kernels available (CUDA, HIP): M=8, N=VECTOR_WIDTH, optional super-clustering (2x2x2 sub-clusters)
- More complex, higher performance on both CPUs and GPUs

### Cluster-Pair Data Structures

**Atom struct** (`src/clusterpair/atom.h`):
```c
typedef struct {
    int natoms;           // Number of atoms in this cluster
    MD_FLOAT bbminx, bbmaxx, bbminy, bbmaxy, bbminz, bbmaxz;  // Bounding box
} Cluster;

typedef struct {
    int nclusters;        // Number of clusters in supercluster
    MD_FLOAT bbminx, ..., bbmaxz;  // Bounding box
} SuperCluster;

typedef struct {
    int Natoms, Nlocal, Nghost, Nmax;
    int Nclusters_local, Nclusters_ghost, Nclusters_max, NmaxGhost, ncj;
    
    // Atomic-level data (AOS or SOA)
    MD_FLOAT *x, *y, *z;  // Positions (separate arrays for SOA, or packed for AOS)
    MD_FLOAT *vx, *vy, *vz;  // Velocities
    int* type;           // Atom type
    
    // Cluster-organized data
    MD_FLOAT* cl_x;      // Positions in cluster format (faster access)
    MD_FLOAT* cl_v;      // Velocities in cluster format
    MD_FLOAT* cl_f;      // Forces in cluster format
    int* cl_t;           // Types in cluster format
    
    // Per-cluster LJ parameters (geometric combination rule optimization)
    MD_FLOAT* cl_sqrt_epsilon;   // sqrt(epsilon) per atom, cluster layout
    MD_FLOAT* cl_sigma3;         // sigma^3 per atom, cluster layout
    
    Cluster* iclusters;         // Local cluster info
    Cluster* jclusters;         // Neighbor cluster info (local + ghost)
    SuperCluster* siclusters;   // Optional super-cluster info
    int* cluster_bin;           // Atom-to-cluster mapping
    MD_UINT* exclusion_filter;  // Masks for intra-cluster exclusions
    
    // Neighbor list and geometry
    int* PBCx, *PBCy, *PBCz;    // Per-atom PBC flags
    MD_FLOAT* sqrt_epsilon_per_type;
    MD_FLOAT* sigma3_per_type;
    Box mybox;                   // Subdomain box
} Atom;
```

**Neighbor struct** (`src/clusterpair/neighbor.h`):
```c
typedef struct {
    int every;                  // Rebuild every N steps
    int ncalls;                 // Number of times built
    int maxneighs;              // Max neighbors per i-cluster
    int* numneigh;              // Number of neighbors per i-cluster
    int* numneigh_masked;       // For SIMD with exclusion masks
    int* numneigh_inner;        // Inner region neighbors
    int half_neigh;             // Half neighbor list flag
    int* neighbors;             // j-cluster indices (AOS or SOA layout)
    unsigned int* neighbors_imask;  // Interaction masks (exclusions/diagonal)
    // MPI fields for shell method
    int Nshell;
    int* numNeighShell;
    int* neighshell;
    int* listshell;
} Neighbor;
```

**Memory layout macros** (cluster-pair):
- Cluster coordinates: `CL_X_INDEX(i)`, `CL_Y_INDEX(i)`, `CL_Z_INDEX(i)` compute offsets in cl_x array
- For 4xN kernel with M >= N: coordinates stored as three contiguous blocks of CLUSTER_M elements
- For 2xNN and GPU kernels: layout is more complex; see force.h for full macro definitions

### Common Parameters

**Runtime parameters** (`src/common/parameter.h`):

Key fields in `Parameter` struct:
- `int ntimes`: Number of timesteps
- `MD_FLOAT dt, dtforce`: Timestep, force timestep (0.5*dt)
- `int nx, ny, nz`: System size (unit cells)
- `MD_FLOAT cutforce, skin, outer_skin`: Force cutoff, skin for neighbor rebuild, additional skin
- `MD_FLOAT cutneigh`: Neighbor list cutoff (cutforce + skin + outer_skin)
- `MD_FLOAT epsilon, sigma, sigma6`: LJ parameters (global or per-type)
- `MD_FLOAT* epsilon_per_type, sigma_per_type`: Per-type LJ parameters (read from types_file)
- `int ntypes`: Number of atom types
- `int LJ_COMB_RULE`: Combination rule (0=single, 1=geometric, 2=full)
- `int half_neigh`: Use half neighbor lists (Newton's 3rd law)
- `int reneigh_every`: Rebuild neighbor list every N steps
- `int nstat`: Print statistics every N steps
- `int super_clustering`: Enable GPU super-clustering

Parameters are loaded from a text file (see `readParameter`, `readTypesFile` in `src/common/parameter.c`). File format is `key value` pairs; types file is two columns (epsilon, sigma).

### LJ Force Computation

**Kernel variants** (`src/clusterpair/force_lj.c`):

- `computeForceLJRef()`: Reference scalar kernel (USE_REFERENCE_KERNEL)
- `computeForceLJ4xnHalfNeigh()`, `computeForceLJ4xnFullNeigh()`: SIMD 4xN kernels
- `computeForceLJ2xnnHalfNeigh()`, `computeForceLJ2xnnFullNeigh()`: SIMD 2xNN kernels
- `computeForceLJCuda()`, `computeForceLJCudaSup()`: CUDA kernels (with/without super-clustering)

**Combination rules** (compile-time via -DLJ_COMB_RULE):
- `LJ_COMB_SINGLE=0`: Single type; broadcast global epsilon/sigma (fastest, no type lookup)
- `LJ_COMB_GEOM=1`: Per-type geometric: epsilon_ij = sqrt(eps_i * eps_j), sigma6_ij = sigma3_i * sigma3_j
- `LJ_COMB_FULL=2`: Full type-pair matrix lookup (not SIMD-optimized)

**Force calculation** (reference kernel):
```c
MD_FLOAT sr2 = 1.0 / rsq;
MD_FLOAT sr6 = sr2 * sr2 * sr2 * sigma6;
MD_FLOAT force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
```
Interaction only if rsq < cutforcesq. Half neighbor-lists apply Newton's 3rd law.

**Verlet list kernel** is much simpler: per-atom iteration with direct neighbor array lookup.

## Testing

### Unit Tests

Compiled into a single binary `mdbench-tests` (not part of default make target):
```bash
make mdbench-tests
./mdbench-tests
```

Unit tests cover:
- `test_parameter.c`: Parameter file parsing
- `test_atom.c`: Atom creation, data layout
- `test_force.c`: Force kernel correctness
- `test_neighbor.c`: Neighbor list building
- `test_integrate.c`: Velocity Verlet integration
- `test_box.c`: Box geometry
- `test_thermo.c`: Thermodynamic properties

Test runner is a minimal header (`tests/test_runner.h`) with `ASSERT_*` macros. Tests are standalone functions returning 0 on pass, 1 on fail.

### Regression Tests (Shell Scripts)

Executed by `make test`:

- `tests/sim_argon_regression.sh`: Run argon simulation, check final T and P
- `tests/sim_copper_fcc_regression.sh`: Run copper FCC simulation
- `tests/regression_energy_lj.sh`: Test energy conservation
- `tests/regression_lj_table.sh`: Compare tabulated vs analytic LJ forces
- `tests/test_double_cutoff.sh`: Test two-tier cutoff (outer_skin)
- `tests/regression_scheme_equiv.sh`: Verify verletlist and clusterpair give same results
- `tests/test_simd_vs_scalar.sh`: Compare SIMD and scalar kernels
- `tests/test_half_neigh.sh`: Verify half-neighbor lists are correct
- `tests/test_data_layout.sh`: Compare AOS vs SOA atom layout
- `tests/test_mpi.sh`: MPI parallel correctness

Most scripts build multiple binary configurations and compare outputs.

### Running Tests

```bash
# Full test suite
make test

# Just unit tests
./mdbench-tests

# Just one regression test
bash tests/sim_argon_regression.sh ./MDBench-VL-GCC-X86-AVX512-DP

# Build and test with specific config
make TOOLCHAIN=GCC OPT_SCHEME=clusterpair SIMD=AVX512 DATA_TYPE=DP test
```

### Test Data

Simulations use input files in `data/`:
- `data/argon/input.gro`: Argon gas configuration (GROMACS format)
- `data/argon/mdbench_params.conf`: Parameter file
- `data/copper_melting/`: Copper FCC system

## Command Quick Reference

```bash
# Configure and build
make TOOLCHAIN=GCC OPT_SCHEME=verletlist SIMD=AVX2 DATA_TYPE=DP

# Clean and rebuild
make clean && make

# Rebuild everything (all configs)
make cleanall

# Code quality
make format                    # Auto-format with clang-format
make asm                       # Generate assembly

# Testing
make test                      # Full test suite
./mdbench-tests                # Unit tests only
make SIMD=NONE test            # Test scalar variant

# Run simulation directly
./MDBench-VL-GCC-X86-AVX512-DP -i data/argon/input.gro -p data/argon/mdbench_params.conf -n 100

# Get compiler info (useful for logging)
make info
```

## Key Files for Development

| File | Purpose |
|------|---------|
| `src/common/parameter.h`, `parameter.c` | Runtime parameter struct and file parsing |
| `src/clusterpair/atom.h` | Cluster-pair atom/cluster data structures |
| `src/verletlist/atom.h` | Verlet list atom structure (simpler) |
| `src/clusterpair/neighbor.h` | Cluster-pair neighbor list structure |
| `src/clusterpair/force.h` | Kernel macro definitions and cluster indexing |
| `src/clusterpair/force_lj.c` | Cluster-pair LJ force kernels (reference, SIMD, CUDA) |
| `src/verletlist/force_lj.c` | Verlet list LJ force kernels |
| `src/common/ljtable.c` | Tabulated LJ force lookup |
| `src/clusterpair/main.c` | Main program entry (argument parsing, simulation loop) |
| `.clang-format` | Code formatting rules (90-char limit, WebKit style) |
| `config.mk` | Build configuration variables |
| `Makefile` | Main build rules (pattern rules, dependency generation, test targets) |

## Important Constraints & Notes

- **Macro-heavy**: cluster-pair code uses extensive macros for cluster indexing and coordinate access (CL_X_INDEX, CI_BASE_INDEX, etc.). These depend on kernel variant, so understand which macros apply to your change.
- **Multiple kernels**: Don't assume a single force kernel; changes may need to be applied to reference, 4xN, 2xNN, CUDA, and GPU super-clustering variants.
- **Precision-aware**: Code switches between float and double via MD_FLOAT. Ensure precision-sensitive operations use correct types.
- **SIMD intrinsic dependencies**: SIMD kernels require specific ISA (e.g., AVX512 kernel won't compile with SIMD=NONE). Test across multiple SIMD levels if adding intrinsic code.
- **GPU code**: CUDA/HIP code paths are orthogonal. NVCC and HIPCC change kernel selection and data layout macros.
- **Per-type parameters**: Geometric combination rule requires per-type sqrt(epsilon) and sigma^3; scalar kernels also need full type-pair epsilon/sigma6 matrices.

