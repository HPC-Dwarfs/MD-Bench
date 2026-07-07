CC = gcc
LINKER = $(CC)

ifeq ($(strip $(ENABLE_MPI)),true)
    CC = mpicc
    DEFINES += -D_MPI
endif

ifeq ($(strip $(ENABLE_OPENMP)),true)
OPENMP = -fopenmp
endif

PROFILE = #-fprofile-generate -fauto-profile -g
ANSI_CFLAGS = #-ansi -pedantic -Wextra

# ARM SIMD options
ifeq ($(strip $(ISA)),ARM)
OPTS = -Ofast
# instead of armvXY-a, we would prefer native, but GCC does not support further
# extensions to it
ifeq ($(strip $(SIMD)),SVE2)
OPTS += -march=armv8.5-a+sve2
endif
ifeq ($(strip $(SIMD)),SVE)
OPTS += -march=armv8.5-a+sve
endif
ifeq ($(strip $(SIMD)),NEON)
OPTS += -march=armv8.5-a
endif
ifeq ($(strip $(SIMD)),NONE)
OPTS += -march=armv8.5-a+nosimd
endif
ASFLAGS =
endif

# X86 SIMD options
ifeq ($(strip $(ISA)),X86)
OPTS = -Ofast -ffast-math -funroll-loops
DEFINES    += -DNO_ZMM_INTRIN
ifeq ($(SIMD),AVX512)
OPTS +=  -march=x86-64-v4 -mavx512f
endif
ifeq ($(strip $(SIMD)),AVX2)
OPTS += -march=x86-64-v3 -mavx2
endif
ifeq ($(strip $(SIMD)),AVX)
OPTS += -march=x86-64-v3 -mno-avx2 -mno-bmi -mno-bmi2 -mno-fma4
endif
ifeq ($(strip $(SIMD)),SSE)
OPTS += -march=x86-64-v2
endif
ifeq ($(strip $(SIMD)),NONE)
OPTS += -fno-tree-loop-vectorize -fno-tree-vectorize -fno-tree-slp-vectorize
endif
ASFLAGS = -masm=intel
endif

CFLAGS = $(PROFILE) $(OPENMP) $(OPTS) -std=c11 $(ANSI_CFLAGS)
LFLAGS = $(PROFILE) $(OPENMP) $(OPTS)
DEFINES += -D_GNU_SOURCE
INCLUDES =
LIBS = -lm
