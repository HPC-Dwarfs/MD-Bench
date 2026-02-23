CC = icx
LINKER = $(CC)

ifeq ($(strip $(ENABLE_MPI)),true)
    CC = mpiicx
    DEFINES += -D_MPI
endif

ifeq ($(strip $(ENABLE_OPENMP)),true)
OPENMP      = -qopenmp
endif

PROFILE     = #-profile-functions -g  -pg

# SIMD options
OPTS        = -O3 -ffast-math
ifeq ($(strip $(SIMD)),AVX512)
OPTS        += -xCORE-AVX512 -qopt-zmm-usage=high
endif
ifeq ($(strip $(SIMD)),AVX2)
OPTS       += -xCORE-AVX2
endif
ifeq ($(strip $(SIMD)),AVX)
OPTS       += -xAVX
endif
ifeq ($(strip $(SIMD)),SSE)
OPTS       += -xSSE4.2
endif
ifeq ($(strip $(SIMD)),NONE)
OPTS        = -O1 -no-vec
endif

DEFINES    += -DNO_ZMM_INTRIN
CFLAGS      = $(PROFILE) $(OPENMP) $(OPTS) -std=c11 #-pedantic-errors
ASFLAGS     = -masm=intel
LFLAGS      = $(PROFILE) $(OPENMP) $(OPTS)
DEFINES    += -D_GNU_SOURCE
INCLUDES    = 
LIBS       += -lm
