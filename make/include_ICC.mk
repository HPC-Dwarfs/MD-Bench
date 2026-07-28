CC = icc
LINKER = $(CC)

ifeq ($(strip $(ENABLE_MPI)),true)
    CC = mpiicc
    DEFINES += -D_MPI
endif

ifeq ($(strip $(ENABLE_OPENMP)),true)
OPENMP      = -qopenmp
else
OPENMP      = -qopenmp-simd
endif

PROFILE     = #-profile-functions -g  -pg

# SIMD options
OPTS        = -Ofast
ifeq ($(strip $(SIMD)),AVX512)
OPTS       += -xCORE-AVX512 -qopt-zmm-usage=high #-debug -g
else
DEFINES    += -DNO_ZMM_INTRIN
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
OPTS       += -no-vec
endif

CFLAGS      = $(PROFILE) $(OPENMP) $(OPTS) -std=c11 -restrict -diag-disable=10441 #-pedantic-errors
ASFLAGS     = -masm=intel
LFLAGS      = $(PROFILE) $(OPENMP) $(OPTS) -diag-disable=10441
DEFINES    += -D_GNU_SOURCE
INCLUDES    =
LIBS        = -lm
