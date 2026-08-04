# Compiler tool chain (GCC/CLANG/ICC/ICX/ONEAPI/NVCC/HIPCC)
TOOLCHAIN ?= GCC
# ISA of instruction code (X86/ARM)
ISA ?= X86
# Instruction set for instrinsic kernels (NONE/<X86-SIMD>/<ARM-SIMD>)
# with X86-SIMD options: NONE/SSE/AVX/AVX_FMA/AVX2/AVX512
# with ARM-SIMD options: NONE/NEON/SVE/SVE2 (SVE not width-agnostic yet!)
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

# Sort atoms at a separate frequency (true or false)
SORT_ATOMS ?= false
# Index variable for tabulated LJ forces (r/rsq)
# r:   uniform grid in distance (matches GROMACS), needs a sqrt per pair
# rsq: uniform grid in squared distance, avoids the sqrt (load-bound variant)
LJ_TABLE_INDEX ?= rsq
# LJ combination rule (single/geometric/full)
# single: single atom type, broadcast global params (fastest, no type lookup)
# geometric: per-type params with geometric combination (default)
# full: full type-pair matrix lookup (not supported in SIMD kernels)
LJ_COMB_RULE ?= geometric
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

#Feature options
OPTIONS =  -DALIGNMENT=64
#OPTIONS +=  More options

################################################################
# DO NOT EDIT BELOW !!!
################################################################
DEFINES =
NBLIST_DATA_LAYOUT_DEFAULT=AOS

ifeq ($(strip $(TOOLCHAIN)),HIPCC)
	VECTOR_WIDTH=1
	SIMD=NONE
	USE_REFERENCE_KERNEL=true
	NBLIST_DATA_LAYOUT_DEFAULT=SOA
endif
ifeq ($(strip $(TOOLCHAIN)),NVCC)
	VECTOR_WIDTH=1
	SIMD=NONE
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

ifeq ($(strip $(SORT_ATOMS)),true)
    DEFINES += -DSORT_ATOMS
endif

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

ifeq ($(strip $(LJ_TABLE_INDEX)),rsq)
    DEFINES += -DLJ_TABLE_RSQ
else ifneq ($(strip $(LJ_TABLE_INDEX)),r)
    $(error Invalid LJ_TABLE_INDEX, must be one of: r, rsq)
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
endif
TAG_SUFFIX = -$(strip $(LJ_COMB_RULE))

ifeq ($(strip $(SIMD)),NONE)
		TOOL_TAG = $(TOOLCHAIN)-$(ISA)
else
		TOOL_TAG = $(TOOLCHAIN)-$(ISA)-$(SIMD)
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
