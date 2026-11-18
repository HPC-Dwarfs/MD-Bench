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
# AOS="atom"-major, SOA="neighbor"-major, CSR=compact (no padding)
# For CPU, auto=AOS; For GPU, auto=SOA
NBLIST_DATA_LAYOUT ?= auto
# Debug
DEBUG ?= false

# Index variable for tabulated LJ forces (r/rsq)
# r:   uniform grid in distance (matches GROMACS), needs a sqrt per pair
# rsq: uniform grid in squared distance, avoids the sqrt (load-bound variant)
LJ_TABLE_INDEX ?= rsq
# LJ combination rule (single/geometric/full)
# single: single atom type, broadcast global params (fastest, no type lookup)
# geometric: per-type params with geometric combination (default)
# full: full type-pair matrix lookup (not supported in SIMD kernels)
LJ_COMB_RULE ?= geometric
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
# Map threadIdx.y to cii and threadIdx.x to cjj (true or false)
# If false, use same thread mapping and reduction instructions as Gromacs
SUPERCLUSTER_INVERSE_THREAD_MAPPING ?= true
# Number of thread-block slices the GPU grid is split into along each
# i-supercluster's j-neighbor list for the super-clustering kernel (grid.y).
# 1 slice means each i-supercluster's neighbors are handled by a single
# block, writing its force directly; >1 splits the neighbor list across
# multiple blocks running concurrently (more parallelism/occupancy) that
# combine their partial forces with atomicAdd (more atomic contention).
# Tune based on system size and GPU occupancy.
SUPERCLUSTER_J_SLICES ?= 4
# Use scalar version (and pray for the compiler to vectorize the code properly)
USE_SCALAR_KERNEL ?= false
# Use reference version (for correction and metrics purposes)
USE_REFERENCE_KERNEL ?= false
# Use SIMD intrinsic kernels for force computation (true or false)
USE_SIMD_KERNEL ?= false
# Compress cutoff-passing neighbors into full SIMD registers before running
# the expensive LJ force stage, instead of masking out failing lanes
# (requires USE_SIMD_KERNEL=true; AVX2/AVX512 only for now)
SIMD_COMPRESS ?= false
# Use SIMD intrinsics to build the Verlet-list neighbor lists (true or false).
# Independent of USE_SIMD_KERNEL, which only affects the force kernel.
# (AVX512 only for now; not implemented for NBLIST_DATA_LAYOUT=CSR)
USE_SIMD_NEIGHBOR ?= false
# Enable XTC output (a GROMACS file format for trajectories)
XTC_OUTPUT ?= false

# Configurations for CUDA
# Use CUDA host memory to optimize transfers
USE_CUDA_HOST_MEMORY ?= false
# Allocate the host-side buffers that are transferred to/from the GPU every
# reneighbor step (atom positions/velocities, per-type LJ params, neighbor
# lists) as pinned (page-locked) memory instead of regular pageable memory.
# Pinned memory gives faster/more predictable cudaMemcpy transfers and, on
# systems with NUMA balancing enabled, avoids stalls from the OS migrating
# pages that go untouched by the CPU between reneighbor steps.
USE_PINNED_MEMORY ?= true

#Feature options
OPTIONS =  -DALIGNMENT=64
#OPTIONS +=  More options

################################################################
# DO NOT EDIT BELOW !!!
################################################################
DEFINES =
NBLIST_DATA_LAYOUT_DEFAULT=AOS

# The architecture follows from SIMD
ifeq ($(origin ISA),command line)
    $(error ISA is derived from SIMD now, use SIMD=X86_SCALAR or ARM_SCALAR for scalar code)
endif
GPU_TOOLCHAIN := $(filter $(strip $(TOOLCHAIN)),NVCC HIPCC)
ifneq ($(GPU_TOOLCHAIN),)
    override SIMD := NONE
endif
ifeq ($(strip $(SIMD)),NONE)
    override SIMD := $(if $(filter aarch64 arm64,$(shell uname -m)),ARM,X86)_SCALAR
endif
X86_SIMD_SET = X86_SCALAR SSE AVX AVX_FMA AVX2 AVX512
ARM_SIMD_SET = ARM_SCALAR NEON SVE SVE2
ifneq ($(filter $(strip $(SIMD)),$(X86_SIMD_SET)),)
    ISA := X86
else ifneq ($(filter $(strip $(SIMD)),$(ARM_SIMD_SET)),)
    ISA := ARM
else
    $(error Invalid SIMD=$(SIMD), must be one of: NONE $(X86_SIMD_SET) $(ARM_SIMD_SET))
endif
SIMD_TAG := $(strip $(SIMD))
# Scalar code is SIMD=NONE for the toolchain files and the defines below
ifneq ($(filter %_SCALAR,$(SIMD)),)
    override SIMD := NONE
endif

ifneq ($(GPU_TOOLCHAIN),)
	VECTOR_WIDTH=1
	USE_REFERENCE_KERNEL=true
	NBLIST_DATA_LAYOUT_DEFAULT=SOA
endif
ifeq ($(strip $(SIMD)),NONE)
	VECTOR_WIDTH=1
	USE_REFERENCE_KERNEL=true
else
ifeq ($(strip $(ISA)),ARM)
    ifeq ($(strip $(SIMD)),NEON)
        __ISA_NEON__=true
        __SIMD_WIDTH_DBL__=2
    else ifeq ($(strip $(SIMD)),SVE)
        __ISA_SVE__=true
		# needs further specification
        __SIMD_WIDTH_DBL__=2
    else ifeq ($(strip $(SIMD)),SVE2)
        __ISA_SVE__=true
        __ISA_SVE2__=true
        # needs further specification
        __SIMD_WIDTH_DBL__=2
    endif
else
# X86
    ifeq ($(strip $(SIMD)),SSE)
        __ISA_SSE__=true
        __SIMD_WIDTH_DBL__=2
    else ifeq ($(strip $(SIMD)),AVX)
        __ISA_AVX__=true
        __SIMD_WIDTH_DBL__=4
    else ifeq ($(strip $(SIMD)),AVX_FMA)
        __ISA_AVX__=true
        __ISA_AVX_FMA__=true
        __SIMD_WIDTH_DBL__=4
    else ifeq ($(strip $(SIMD)),AVX2)
        __ISA_AVX2__=true
        __SIMD_WIDTH_DBL__=4
    else ifeq ($(strip $(SIMD)),AVX512)
        __ISA_AVX512__=true
        __SIMD_WIDTH_DBL__=8
    endif
endif

# SIMD width is specified in double-precision, hence it needs
# to be adjusted for single-precision cases
ifeq ($(strip $(DATA_TYPE)), SP)
    VECTOR_WIDTH=$(shell echo $$(( $(__SIMD_WIDTH_DBL__) * 2 )))
else
    VECTOR_WIDTH=$(__SIMD_WIDTH_DBL__)
endif
endif
ifeq ($(strip $(ATOM_DATA_LAYOUT)),AOS)
    DEFINES +=  -DATOM_POSITION_AOS
endif
ifeq ($(strip $(NBLIST_DATA_LAYOUT)),auto)
    NBLIST_DATA_LAYOUT=$(NBLIST_DATA_LAYOUT_DEFAULT)
endif
ifeq ($(strip $(NBLIST_DATA_LAYOUT)),AOS)
    DEFINES +=  -DNBLIST_AOS
else ifeq ($(strip $(NBLIST_DATA_LAYOUT)),CSR)
    DEFINES +=  -DNBLIST_CSR
else ifeq ($(strip $(NBLIST_DATA_LAYOUT)),SOA)
    DEFINES +=  -DNBLIST_SOA
else
    $(error Invalid NBLIST_DATA_LAYOUT: $(NBLIST_DATA_LAYOUT). Must be one of: auto, AOS, SOA, CSR)
endif
ifeq ($(strip $(SUPERCLUSTER_DATA_LAYOUT)),AOS3)
    DEFINES +=  -DPOSITION_AOS3_SUP
else ifeq ($(strip $(SUPERCLUSTER_DATA_LAYOUT)),AOS4)
    DEFINES +=  -DPOSITION_AOS4_SUP
else
    DEFINES +=  -DPOSITION_SOA_SUP
endif
ifeq ($(strip $(DATA_TYPE)),SP)
    DEFINES +=  -DPRECISION=1
else
    DEFINES +=  -DPRECISION=2
endif

# Former build options that are runtime options now
RUNTIME_OPTIONS = SORT_ATOMS:--sort-every ACCURATE_RECIPROCAL:--accurate-recip KERNEL_PREFETCH:--no-prefetch
ifeq ($(origin EXCLUSIONS),command line)
    $(error EXCLUSIONS is gone, the molecular exclusions are on when the system has molecules (-gmxbench))
endif
$(foreach o,$(RUNTIME_OPTIONS),$(if $(filter-out undefined,$(origin $(word 1,$(subst :, ,$(o))))),$(error $(word 1,$(subst :, ,$(o))) is a runtime option now, use $(word 2,$(subst :, ,$(o))) of the binary)))

# Translate LJ_COMB_RULE to compiler define
ifeq ($(strip $(LJ_COMB_RULE)),single)
    DEFINES += -DLJ_COMB_RULE=0
else ifeq ($(strip $(LJ_COMB_RULE)),geometric)
    DEFINES += -DLJ_COMB_RULE=1
else ifeq ($(strip $(LJ_COMB_RULE)),full)
    DEFINES += -DLJ_COMB_RULE=2
else
    $(error Invalid LJ_COMB_RULE, must be one of: single, geometric, full)
endif

# Translate COULOMB_RULE to compiler defines
ifeq ($(strip $(COULOMB_RULE)),none)
    DEFINES += -DCOULOMB_RULE=0
else ifeq ($(strip $(COULOMB_RULE)),rf)
    DEFINES += -DCOULOMB_RULE=1
else ifeq ($(strip $(COULOMB_RULE)),ewald)
    DEFINES += -DCOULOMB_RULE=2
else ifeq ($(strip $(COULOMB_RULE)),ewald_table)
    ifneq ($(filter $(strip $(TOOLCHAIN)),NVCC HIPCC),)
        $(error COULOMB_RULE=ewald_table is only supported by the CPU kernels)
    endif
    DEFINES += -DCOULOMB_RULE=3
else
    $(error Invalid COULOMB_RULE, must be one of: none, rf, ewald, ewald_table)
endif

ifeq ($(strip $(LJ_TABLE_INDEX)),rsq)
    DEFINES += -DLJ_TABLE_RSQ
else ifneq ($(strip $(LJ_TABLE_INDEX)),r)
    $(error Invalid LJ_TABLE_INDEX, must be one of: r, rsq)
endif

ifeq ($(strip $(ASM_SYNTAX)),intel)
    ASM_SYNTAX_FLAG = -masm=intel
else ifeq ($(strip $(ASM_SYNTAX)),att)
    ASM_SYNTAX_FLAG = -masm=att
else
    $(error Invalid ASM_SYNTAX, must be one of: intel, att)
endif

ifeq ($(strip $(MEM_TRACER)),true)
    DEFINES += -DMEM_TRACER
endif

ifeq ($(strip $(INDEX_TRACER)),true)
    DEFINES += -DINDEX_TRACER
endif

ifeq ($(strip $(COMPUTE_STATS)),true)
    DEFINES += -DCOMPUTE_STATS
endif

ifeq ($(strip $(XTC_OUTPUT)),true)
    DEFINES += -DXTC_OUTPUT
endif

ifeq ($(strip $(USE_PINNED_MEMORY)),true)
    DEFINES += -DUSE_PINNED_MEMORY
endif

ifeq ($(strip $(USE_SCALAR_KERNEL)),true)
    DEFINES += -DUSE_SCALAR_KERNEL
endif

ifeq ($(strip $(USE_REFERENCE_KERNEL)),true)
    DEFINES += -DUSE_REFERENCE_KERNEL
endif

ifeq ($(strip $(DEBUG)),true)
    DEFINES += -DDEBUG
endif

ifneq ($(VECTOR_WIDTH),)
    DEFINES += -DVECTOR_WIDTH=$(VECTOR_WIDTH)
endif

ifeq ($(strip $(USE_SIMD_KERNEL)),true)
    DEFINES += -D__SIMD_KERNEL__
    ifeq ($(strip $(SIMD_COMPRESS)),true)
        DEFINES += -D__SIMD_COMPRESS__
    endif
    # SVE/SVE2 are vector-length-agnostic; use the single-loop VLA kernel there instead of the fixed-width main+tail split.
    ifneq ($(filter $(strip $(ISA))-$(strip $(SIMD)),ARM-SVE ARM-SVE2),)
        DEFINES += -D__SIMD_VLA__
    endif
endif

ifeq ($(strip $(USE_SIMD_NEIGHBOR)),true)
    DEFINES += -D__SIMD_NEIGHBOR__
endif

ifeq ($(strip $(__SSE__)),true)
    DEFINES += -D__ISA_SSE__
endif

ifeq ($(strip $(__ISA_AVX__)),true)
    DEFINES += -D__ISA_AVX__
endif

ifeq ($(strip $(__ISA_AVX_FMA__)),true)
    DEFINES += -D__ISA_AVX_FMA__
endif

ifeq ($(strip $(__ISA_AVX2__)),true)
    DEFINES += -D__ISA_AVX2__
endif

ifeq ($(strip $(__ISA_AVX512__)),true)
    DEFINES += -D__ISA_AVX512__
endif

ifeq ($(strip $(__ISA_NEON__)),true)
    DEFINES += -D__ISA_NEON__
endif

ifeq ($(strip $(__ISA_SVE__)),true)
    DEFINES += -D__ISA_SVE__
endif

ifeq ($(strip $(__ISA_SVE2__)),true)
    DEFINES += -D__ISA_SVE2__
endif

ifeq ($(strip $(OPT_SCHEME)),clusterpair)
    DEFINES += -DCLUSTER_PAIR
endif

ifeq ($(strip $(SUPERCLUSTER_INVERSE_THREAD_MAPPING)),true)
    DEFINES += -DSUPERCLUSTER_INVERSE_THREAD_MAPPING
endif

DEFINES += -DGPU_J_SLICES=$(strip $(SUPERCLUSTER_J_SLICES))

ifeq ($(strip $(OPT_SCHEME)),verletlist)
		OPT_TAG = VL
    ifeq ($(strip $(USE_SIMD_KERNEL)),true)
        OPT_TAG := $(OPT_TAG)-SIMD
    endif
    ifeq ($(strip $(USE_SIMD_NEIGHBOR)),true)
        OPT_TAG := $(OPT_TAG)-NBSIMD
    endif
else ifeq ($(strip $(OPT_SCHEME)),clusterpair)
    OPT_TAG = CP-$(CLUSTER_PAIR_KERNEL)
    # The C reference kernel of the CPU (forced without SIMD)
    ifeq ($(strip $(USE_REFERENCE_KERNEL))-$(filter $(strip $(TOOLCHAIN)),NVCC HIPCC),true-)
        OPT_TAG = CP-ref
    endif
endif
TAG_SUFFIX = -$(strip $(LJ_COMB_RULE))
ifneq ($(strip $(COULOMB_RULE)),none)
    TAG_SUFFIX := $(TAG_SUFFIX)-$(strip $(COULOMB_RULE))
endif
ifeq ($(strip $(ENABLE_OPENMP)),true)
    TAG_SUFFIX := $(TAG_SUFFIX)-omp
endif
ifeq ($(strip $(ENABLE_MPI)),true)
    TAG_SUFFIX := $(TAG_SUFFIX)-mpi
endif

# GPU builds run no CPU SIMD code
ifneq ($(GPU_TOOLCHAIN),)
    TOOL_TAG = $(TOOLCHAIN)
else
    TOOL_TAG = $(TOOLCHAIN)-$(SIMD_TAG)
endif

ifeq ($(strip $(SIMD_COMPRESS)),true)
		TOOL_TAG := $(TOOL_TAG)-CMP
endif

ifeq ($(strip $(OPT_SCHEME)),clusterpair)
    ifeq ($(strip $(CLUSTER_PAIR_KERNEL)),auto)
        DEFINES += -DCLUSTERPAIR_KERNEL_AUTO
    else ifeq ($(strip $(CLUSTER_PAIR_KERNEL)),supercluster)
        DEFINES += -DCLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
    else ifeq ($(strip $(CLUSTER_PAIR_KERNEL)),gpusimple)
        DEFINES += -DCLUSTERPAIR_KERNEL_GPU_SIMPLE
    else ifeq ($(strip $(CLUSTER_PAIR_KERNEL)),4xN)
        DEFINES += -DCLUSTERPAIR_KERNEL_4XN
    else ifeq ($(strip $(CLUSTER_PAIR_KERNEL)),2xN)
        DEFINES += -DCLUSTERPAIR_KERNEL_2XN
    else ifeq ($(strip $(CLUSTER_PAIR_KERNEL)),2xNN)
        DEFINES += -DCLUSTERPAIR_KERNEL_2XNN
    else
        $(error Invalid CLUSTER_PAIR_KERNEL, must be one of: auto, 4xN, 2xNN, 2xN, gpusimple, supercluster)
    endif
endif
