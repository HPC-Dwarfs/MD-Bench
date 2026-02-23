CC  = nvcc
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

#
# A100 + Native
CFLAGS   = -O3 -arch=sm_80 -march=native -ffast-math -funroll-loops --forward-unknown-to-host-compiler # -fopenmp
# A40 + Native
#CFLAGS   = -O3 -arch=sm_86 -march=native -ffast-math -funroll-loops --forward-unknown-to-host-compiler # -fopenmp
# Cascade Lake
#CFLAGS   = -O3 -march=cascadelake  -ffast-math -funroll-loops --forward-unknown-to-host-compiler # -fopenmp
# For GROMACS kernels, we need at least sm_61 due to atomicAdd with doubles
# TODO: Check if this is required for full neighbor-lists and just compile kernel for that case if not
#CFLAGS   = -O3 -g -arch=sm_61 # -fopenmp
ASFLAGS  =  -masm=intel
LFLAGS   =
DEFINES  += -D_GNU_SOURCE -DCUDA_TARGET=0 -DNO_ZMM_INTRIN  #-DLIKWID_PERFMON
INCLUDES = $(MPI_HOME) $(LIKWID_INC)
LIBS     = -lm -lcuda -lcudart $(LIKWID_LIB) $(MPI_LIB)#-llikwid
