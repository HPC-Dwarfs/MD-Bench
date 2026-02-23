# MD-Bench

| :exclamation:  The repository was restructured and cleaned up and in this process the history was rewritten. Please make fresh forks and clones.   |
|-----------------------------------------|

MD-Bench is a toolbox for the performance engineering of short-range force
calculation kernels on molecular-dynamics applications. It aims at covering
state-of-the-art algorithms from different community codes such as
LAMMPS and GROMACS.

## Getting started

Clone the repository from GitHub:

```shell=
git clone https://github.com/RRZE-HPC/MD-Bench.git
```

Edit config.mk and configure the compiler toolchain

```makefile=
# Compiler tool chain (GCC/CLANG/ICC/ICX/ONEAPI/NVCC)
TOOLCHAIN ?= CLANG
```

Best supported are ICC (deprecated legacy Intel compiler) and ICX (LLVM based
Intel compiler). Choose NVCC to enable CUDA GPU kernels.

The toolchain settings are located in the `./make` directory. Review the
settings for the configured toolchain. You can configure different settings in
`config.mk`, for starters on a X86 based system the defaults are fine.

To build the binary call, don't forget to load the compiler module on the
NHR@FAU clusters (e.g. `module load intel`):

```shell=
make
```

While the Makefile works with any version of GNU make, some features require GNU
make > v4.

You can run MD-Bench without any arguments:

```shell=
./MDBench-VL-ICX-X86-AVX2-DP
```

## Build system

MD-Bench uses a Makefile with pattern rules and automatic dependency generation.
If you add source files you do not need to change the Makefile as long as the
sources are placed either in the `./src/verletlist/`, `./src/clusterpair/` or
`./src/common` directories. If you change a file, all object files that depend
on it are rebuild.

All configuration variables can be overwritten from the command line, e.g. to
build with ICC without changing `./config.mk` build with:

```shell=
make TOOLCHAIN=ICC
```

Multiple configurations can be build at the same time. Every configuration has a
unique binary name `./MDBENCH-<build tag>`. Intermediate build results are
located in a `./build/build-<build tag>/` directory.

All make targets act on the current configuration set in `./config.mk`, but this
can be of course overwritten on the command line.

Supported make targets:

- `make`: Build the binary for current configuration.
- `make clean`: Remove intermediate build results.
- `make distclean`: Remove intermediate build results and binary. Also removes
generated tags and clangd files, more on that later.
- `make cleanall`: Remove all generated files. **Note**: This target applies to
all configurations.
- `make info` Output compiler version, useful for logging in automated benchmark
scripts.
- `make asm`: Generate assembly output of all source files. The assembly files
are placed in the intermediate build directory.
- `make format`: Reformat all source files with `clang-format` using the format
specification in `.clang-format`.

### Build time options

- `TOOLCHAIN`: Determines which toolchain makefile is included
- `ISA`: No usage apart from tag strings
- `SIMD`: Controls the generation of intrinsic kernels for clusterpair
- `OPT_SCHEME`: Algorithmic variant (verletlist or clusterpair), different
source directories and main routines are used
- `ENABLE_LIKWID`: Turn on LIKWID instrumentation, the LIKWID library has to be available
- `DATA_TYPE`: Switch between single precision and double precision floating
point. This is controlled by defines.
- `ATOM_DATA_LAYOUT`: Switch between array-of-structure (AOS) and structure-of-array
(SOA) layout for atom positions and forces. Tradeoff between better cache
utilisation and easier SIMD vectorization.
- `DEBUG`: Enable additional debug output
- `SORT_ATOMS`: Resort atoms to ensure that atoms that are nearby are also close
to each other in the data structures
- `ONE_ATOM_TYPE`: Simulate only one atom type and do not perform table lookup for parameters.
- `ENABLE_OMP_SIMD`: This enforces the use of `#pragma omp simd` for the
verletlist half-neighbour list force kernel. Without is the Intel compiler (at
least ICC) refuses to do SIMD vectorization.
- `USE_REFERENCE_VERSION`: Enforce usage of C implementation for clusterpair
algorithm for validation
- `USE_CUDA_HOST_MEMORY`: Enable pinned host memory for faster host-device transfers
- `ENABLE_MPI:` Turn on the MPI parallel version of the code

### Build for GPU targets

MD-Bench currently only supports Nvidia GPUs using CUDA kernels. To enable CUDA
kernels you need to specify `NVCC` as toolchain. The CUDA source code is in the
same source directories with Cuda suffix and `.cu` as file type ending. If
`NVCC` is set as toolchain, all supported kernels are automatically set to their
CUDA variants at build time. This means a binary either supports CPU kernels or
GPU kernels.

## Command line arguments

MD-Bench can be executed without any arguments, in this case the full neighbor
list testcase with LJ force will be computed for 200 steps and a size of
32x32x32 unit cells.

- `-p / --params <string>`: file to read parameters from (can be specified more
than once). Default initialization sets parameters for default LJ testcase.
*`-f <string>`: force field (lj, eam), default lj. For anything different than
lj you also need to provide spcific parameter file.
- `-i <string>`:  input file with atom positions (dump). MD-Bench supports
Brookhaven protein data bank (.pdb), GROMACS GROMOS87 (.gro), and LAMMPS dump
(.dmp) file formats
- `-e <string>`:  input file for EAM parameters
- `-n / --nsteps <int>`:  set number of timesteps for simulation (default 200)
- `-nx/-ny/-nz <int>`:  set linear dimension of systembox in x/y/z direction
(default 32 in every dimension)
- `-half <int>`: use half (1) or full (0) neighbor lists (default 0 - full
neighbor list)
- `-r / --radius <real>`:   set cutoff radius (default 2.5)
- `-s / --skin <real>`:   set skin (verlet buffer, default 0.3)
- `-w <file>`:  write input atoms to file
- `--freq <real>`:  processor frequency (GHz), used to calculate cycle metrics
(default 2.4)
- `--vtk <string>`:    VTK output file for visualization

## Available testcases

For all variants you can switch between single precision and double precision
and between AOS versus SOA data layouts using build time options. You can use
the half neighbour list algorithm instead of the default full neighbour list by
setting `-half 1`. To enforce SIMD vectorization for the half neighbour list
algorithm you can set the option `ENABLE_OMP_SIMD=true`.

### Lennard-Jones potential for solid copper

Just start without any command line argument, this is the default testcase. You
may change the number of timesteps using the `-n` options and change the problem
size using the `-nz, -ny, -nz` options.

### EAM potential for solid copper

Call MD-Bench as follows:

```shell=
./MDBench-<TAG> -n 400 -f eam -e ./data/Cu_u3.eam
```

Two different EAM variants are available: `Cu_u3.eam` and `Cu_u6.eam`. The EAM
potential is currently only available for verletlist.

### Lennard-Jones potential for melted copper

The melted copper testcase has only 32000 atoms in the default configuration.
Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -n 400  -i ./data/copper_melting/input_lj_cu_one_atomtype_20x20x20.dmp
```

### Lennard-Jones potential for melted copper with explicit types

Compile MD-Bench with `EXPLICIT_TYPES=true` in `config.mk`.

Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -n 400  -i ./data/copper_melting/input_lj_cu_one_atomtype_20x20x20.dmp
```

**This testcase currently segvaults!**

### EAM potential for melted copper

Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -n 400 -f eam -e ./data/Cu_u3.eam  -i ./data/copper_melting/input_eam_cu_one_atomtype_20x20x20.dmp
```

Two different EAM variants are available: `Cu_u3.eam` and `Cu_u6.eam`. The EAM
potential is currently only available for verletlist.t.

### Lennard-Jones potential for argon gas

Call MD-Bench as follows:

```shell=bash
./MDBench-<TAG> -i ./data/argon/input.gro  -p ./data/argon/mdbench_params.conf
```

## Citations

Rafael Ravedutti Lucio Machado, Jan Eitzinger, Jan Laukemann, Georg Hager,
Harald Köstler and Gerhard Wellein: MD-Bench: A performance-focused prototyping
harness for state-of-the-art short-range molecular dynamics algorithms. Future
Generation Computer Systems
([FGCS](https://www.sciencedirect.com/journal/future-generation-computer-systems)),
Volume 149, 2023, Pages 25-38, ISSN 0167-739X, DOI:
[https://doi.org/10.1016/j.future.2023.06.023](https://doi.org/10.1016/j.future.2023.06.023)

Rafael Ravedutti Lucio Machado, Jan Eitzinger, Harald Köstler, and Gerhard
Wellein: MD-Bench: A generic proxy-app toolbox for state-of-the-art molecular
dynamics algorithms. Accepted for [PPAM](https://ppam.edu.pl/) 2022, the 14th
International Conference on Parallel Processing and Applied Mathematics, Gdansk,
Poland, September 11-14, 2022. PPAM 2022 Best Paper Award. Preprint:
[arXiv:2207.13094](https://arxiv.org/abs/2207.13094), DOI:
[https://dl.acm.org/doi/10.1007/978-3-031-30442-2_24](https://dl.acm.org/doi/10.1007/978-3-031-30442-2_24)

## Credits

MD-Bench is developed by the Erlangen National High Performance Computing Center
([NHR@FAU](https://hpc.fau.de/)) at the University of Erlangen-Nürnberg.

## License

[LGPL-3.0](https://github.com/RRZE-HPC/MD-Bench/blob/master/LICENSE)
