# Compiler tool chain (GCC/CLANG/ICC/ICX/ONEAPI/NVCC/HIPCC)
TOOLCHAIN ?= GCC
# Assembly syntax used when emitting .s files on X86 (intel/att), passed to the
# compiler as -masm=<value>. Only affects the X86 toolchains (GCC/ICC/ICX/CLANG/
# NVCC/HIPCC); ARM's assembler has no equivalent flag, so this is ignored there.
ASM_SYNTAX ?= intel
# Instruction set, it also selects the architecture (ISA) of the build:
# X86: X86_SCALAR/SSE/AVX/AVX_FMA/AVX2/AVX512
# ARM: ARM_SCALAR/NEON/SVE/SVE2 (SVE not width-agnostic yet!)
# *_SCALAR is plain scalar code without SIMD extensions, NONE the one of the build machine
SIMD ?= AVX512
# Optimization scheme (verletlist/clusterpair)
OPT_SCHEME ?= verletlist
# Enable likwid (true or false)
ENABLE_LIKWID ?= false
# Enable OpenMP parallelization (true or false)
ENABLE_OPENMP ?= true
# Enable MPI parallelization
ENABLE_MPI ?= false
# SP or DP
DATA_TYPE ?= SP
# AOS or SOA
ATOM_DATA_LAYOUT ?= AOS
# Neighbor-lists data layout (auto/AOS/SOA/CSR)
# AOS="atom"-major, SOA="neighbor"-major, CSR=compact (no padding, CPU only)
# For CPU, auto=AOS; for GPU (NVCC, HIPCC), auto=AOS with clusterpair and SOA with VL
NBLIST_DATA_LAYOUT ?= auto
# Debug
DEBUG ?= false

# Electrostatics (none/rf/ewald/ewald_table), plain LJ kernels of both schemes
# none: LJ only (default)
# rf: reaction-field with infinite dielectric, shifted at the cut-off
# ewald: real-space Ewald part with an analytical force correction
# ewald_table: real-space Ewald part with GROMACS' tabulated correction (CPU kernels)
COULOMB_RULE ?= none
# Trace memory addresses for cache simulator (true or false)
MEM_TRACER ?= false
# Trace indexes and distances for gather-md (true or false)
INDEX_TRACER ?= false
# Compute statistics
COMPUTE_STATS ?= false

# Configurations for clusterpair optimization scheme
# Cluster pair kernel variant (auto/4xN/2xNN/2xN/gpusimple/supercluster)
CLUSTER_PAIR_KERNEL ?= auto
# Data layout for super-clustering kernels (AOS3/AOS4/SOA)
SUPERCLUSTER_DATA_LAYOUT ?= AOS3
# Use scalar version (and pray for the compiler to vectorize the code properly)
USE_SCALAR_KERNEL ?= false
# Use reference version (for correction and metrics purposes)
USE_REFERENCE_KERNEL ?= false
# Enable XTC output (a GROMACS file format for trajectories)
XTC_OUTPUT ?= false

# Clamp the squared distance of the SIMD cluster-pair kernels to MIN_DISTANCE_SQ after
# the cut-off check, as GROMACS does. The padding slots are masked by their atom id, so
# the clamp only matters for real atoms closer than sqrt(MIN_DISTANCE_SQ).
CLUSTER_RSQ_CLAMP ?= true

#Feature options
OPTIONS =  -DALIGNMENT=64
#OPTIONS +=  More options
