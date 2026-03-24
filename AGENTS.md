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

# Build with CUDA support
make TOOLCHAIN=NVCC

# Clean build artifacts
make clean

# Clean everything including binaries
make cleanall

# Format source code
make format

# Generate assembly output
make asm

# Show compiler information
make info
```

### Configuration Options

Key configuration options in `config.mk`:

- `TOOLCHAIN`: Compiler toolchain (GCC/CLANG/ICC/ICX/NVCC/HIPCC)
- `OPT_SCHEME`: Algorithmic variant (verletlist/clusterpair)
- `SIMD`: SIMD instruction set (NONE/SSE/AVX/AVX2/AVX512/NEON/SVE/SVE2)
- `DATA_TYPE`: Precision (SP for single, DP for double)
- `ATOM_DATA_LAYOUT`: Data layout (AOS for array-of-structures, SOA for structure-of-arrays)
- `ENABLE_MPI`: Enable MPI parallelization
- `ENABLE_OPENMP`: Enable OpenMP parallelization

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

### Kernel Variants

For cluster-pair implementation, several kernel variants are available:
- Reference kernel (scalar)
- SIMD 4xN kernel 
- SIMD 2xNN kernel
- CUDA GPU kernel
- CUDA super-cluster kernel

The kernel selection is controlled through build-time options and affects the generated binary name.

## Testing and Execution

### Running the Application

Execute with default parameters:
```bash
./MDBench-<TAG>
```

Common command line arguments:
- `-n <steps>`: Number of timesteps (default 200)
- `-nx/-ny/-nz <int>`: System dimensions (default 32x32x32)
- `-f <field>`: Force field (lj for Lennard-Jones, eam for EAM)
- `-half <0|1>`: Use half (1) or full (0) neighbor lists
- `-r <radius>`: Cutoff radius (default 2.5)
- `-s <skin>`: Skin for Verlet buffer (default 0.3)

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

## Important Gotchas

1. **Build Tags**: Each configuration produces a uniquely named binary with a tag reflecting the build options. Multiple configurations can coexist.

2. **CUDA Limitations**: CUDA builds only support GPU kernels; CPU kernels are disabled when TOOLCHAIN=NVCC.

3. **SIMD Dependencies**: SIMD kernels require specific hardware support. The reference kernel is used automatically when SIMD is not available.

4. **Data Layout Impact**: AOS vs SOA layout affects both performance and memory access patterns significantly.

5. **Neighbor List Construction**: The choice between full and half neighbor lists affects both correctness and performance.

6. **Precision Effects**: Single vs double precision not only affects accuracy but also vectorization efficiency.

7. **Header Dependencies**: The force.h header contains complex macro definitions that control kernel generation.

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
- Struct typedefs end with `_t` suffix
- File-local functions use `static` keyword
- Extensive use of preprocessor conditionals for feature selection
- Heavy use of macros for code generation in performance-critical kernels

The project uses `.clang-format` for code formatting consistency.