CC  = hipcc
LINKER = $(CC)

ifeq ($(strip $(ENABLE_MPI)),true)
DEFINES += -D_MPI
MPI_LIB  = -L/apps/SPACK/0.19.1/opt/linux-almalinux8-zen/nvhpc-23.7/openmpi-4.1.6-ojusv6lrh7e5o7ktibh2qaj2yuzxyzeg/lib \
           -L/apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/hwloc-2.8.0-bneg6wh22jt37qyr2hghz5vmrdk6txyt/lib \
           -L/apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/ucx-1.13.1-woaymodwh7p66njpgt76d7fyqyv7srl3/lib \
           -L/apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/libevent-2.1.12-rzghn5u3lsuysqbbf7xqwq4kg4gkr7mf/lib \
           -L/apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/pmix-4.1.3-lksytzxfs4yylwzgc7dlcihvkqnc7e4d/lib \
           -L/usr/lib64 \
           -Xlinker -rpath -Xlinker /apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/nvhpc-23.7-bzxcokzjvx4stynglo4u2ffpljajzlam/Linux_x86_64/23.7/compilers/lib \
           -Xlinker -rpath -Xlinker /apps/SPACK/0.19.1/opt/linux-almalinux8-zen/nvhpc-23.7/openmpi-4.1.6-ojusv6lrh7e5o7ktibh2qaj2yuzxyzeg/lib \
           -Xlinker -rpath -Xlinker /apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/hwloc-2.8.0-bneg6wh22jt37qyr2hghz5vmrdk6txyt/lib \
           -Xlinker -rpath -Xlinker /apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/ucx-1.13.1-woaymodwh7p66njpgt76d7fyqyv7srl3/lib \
           -Xlinker -rpath -Xlinker /apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/libevent-2.1.12-rzghn5u3lsuysqbbf7xqwq4kg4gkr7mf/lib \
           -Xlinker -rpath -Xlinker /apps/SPACK/0.19.1/opt/linux-almalinux8-zen/gcc-8.5.0/pmix-4.1.3-lksytzxfs4yylwzgc7dlcihvkqnc7e4d/lib \
           -lmpi
MPI_HOME = -I/apps/SPACK/0.19.1/opt/linux-almalinux8-zen/nvhpc-23.7/openmpi-4.1.6-ojusv6lrh7e5o7ktibh2qaj2yuzxyzeg/include
endif

ANSI_CFLAGS  = -ansi
ANSI_CFLAGS += -std=c99
ANSI_CFLAGS += -pedantic
ANSI_CFLAGS += -Wextra

# Target GPU architecture must be provided explicitly on the command line or in the
# environment, e.g.: make GPU_ARCH=gfx90a
# Common values: gfx90a = MI250X, gfx940 = MI300A (APU), gfx942 = MI300X
ifeq ($(strip $(GPU_ARCH)),)
$(error GPU_ARCH is not set. Please specify a supported AMD GPU target, e.g. 'make GPU_ARCH=gfx90a')
endif

# CPU-side routines that never touch the GPU (e.g. buildNeighbor()) are
# still OpenMP-parallelized in source (#pragma omp, guarded by _OPENMP);
# forward -fopenmp to the host compiler so those pragmas actually take
# effect instead of silently compiling as serial code.
ifeq ($(strip $(ENABLE_OPENMP)),true)
OPENMP = -fopenmp
endif

CFLAGS   = -O3 --offload-arch=$(GPU_ARCH) -march=native -ffast-math -funroll-loops $(OPENMP)

ASFLAGS  =  -masm=intel
LFLAGS   = $(OPENMP)
INCLUDES = $(LIKWID_INC) $(MPI_HOME) -I/opt/rocm/include
DEFINES  += -D_GNU_SOURCE -DCUDA_TARGET=1 -DNO_ZMM_INTRIN #-DLIKWID_PERFMON
LIBS     = -lm $(LIKWID_LIB) $(MPI_LIB) -lamdhip64 -lroctx64 #-llikwid
