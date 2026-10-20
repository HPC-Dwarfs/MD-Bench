# MD-Bench vs. the GROMACS non-bonded benchmark

`gmx nonbonded-benchmark` (`src/gromacs/nbnxm/benchmark` in GROMACS) times the NBNxM
pair kernels on 1000 SPC/E water molecules (3000 atoms, replicated `-size` times) with a
single fixed pair list. MD-Bench runs the same system and the same kernel cluster
geometry (4xM with M = SIMD width) through its generic cluster-pair 4xN kernels, with
no kernel specific to the benchmark:

```shell
make OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL=4xN COULOMB_RULE=ewald ACCURATE_RECIPROCAL=true
OMP_NUM_THREADS=1 ./MDBench-CP-4xN-<...>-ewald-accrecip -gmxbench -size 1 -half 1 --kernel-bench 500
```

| Feature | GROMACS benchmark | MD-Bench |
|---|---|---|
| system | `bench_coords.h`, put in the box atom by atom | `src/common/gmxcoords.h` (same numbers), same wrapping |
| cut-off, list buffer | 1.0 nm, no buffer | `-gmxbench` defaults to `-r 1.0 -s 0` |
| Coulomb | `-coulomb ewald` (analytical correction) / `reaction-field` | `COULOMB_RULE=ewald` / `rf` |
| exclusions | whole molecule, pair list masks | molecule index per atom, masked part of the list |
| LJ | geometric rule, no half-LJ | `LJ_COMB_RULE=geometric` |
| result | timing only | `--kernel-bench N` (time per call, pairs/us) |

## Files

- `nbnxm_reference.py`: exact float64 O(N^2) forces of the system (erf/erfc, minimum
  image), the oracle that neither code shares approximations with.
- `compare_forces.py`: error of a force dump relative to the RMS reference force.
- `gromacs-force-dump.patch`: adds `gmx nonbonded-benchmark -forces <file>` (the
  benchmark prints only timings); apply it to a GROMACS tree with `git apply`.
- `run_comparison.sh`: forces and kernel times, MD-Bench vs. GROMACS.
- `../test_coulomb_exclusions.sh`: MD-Bench against the reference, part of `make test`.

## Running the comparison

```shell
git -C ~/gromacs apply ~/MD-Bench/tests/gromacs_nonbonded/gromacs-force-dump.patch
# build gmx (cmake -DGMX_SIMD=AVX2_256 -DGMX_GPU=OFF -DGMX_MPI=OFF -DGMX_BUILD_OWN_FFTW=ON ...)
GMX=~/gromacs/build/bin/gmx PIN_CORE=3 tests/gromacs_nonbonded/run_comparison.sh 1 8
```

## Things to know

- **The benchmark's LJ is C6/6 and C12/12.** The benchmark hands `c6Oxygen`/`c12Oxygen`
  to the NBNxM atom data, which expects the force field's pre-multiplied `6*C6` and
  `12*C12`. The force it computes is `(C12 r^-12 - C6 r^-6) / r^2`, a LJ potential with
  C6/6 and C12/12 (much weaker LJ forces than real SPC/E). MD-Bench and the
  reference use these effective parameters so the forces match.
- **`ACCURATE_RECIPROCAL=true` is needed for matching forces.** MD-Bench's LJ kernels
  use the hardware reciprocal estimate (12-14 bits, also in DP), which gives 1e-3
  relative force errors. GROMACS iterates its estimate to full precision.
- **Half lists need one thread.** The half-list kernels scatter j-forces with
  plain stores from all threads, so with several OpenMP threads the forces are wrong and
  vary from run to run (4 threads, 3000 atoms: RMS error 5e-2 between two runs). Full
  lists are unaffected. This is how the existing 4xN half kernel behaves.
- **Pair-list quality differs.** GROMACS' grid gives 1.34M atom pairs for the 3000-atom
  system, MD-Bench's binning 2.03M. The kernel throughput in computed pairs per us is
  therefore the comparable number; time per call and useful pairs/us also include the
  list quality.
- **`-gmxbench` is a kernel benchmark, not a stable MD system.** The molecules have no
  bonded or constraint terms and the hydrogens have no LJ, so integrating it for many
  steps collapses (the trajectory is nevertheless independent of the neighbor-list
  cadence and list type, which exercises the exclusion handling across re-neighboring).
- Supported SIMD back-ends for exclusions: AVX2, AVX512, NEON, SVE (only AVX2 and AVX512
  were built, only AVX2 was run); plain AVX and the scalar emulation lack the integer mask
  operations.
- Not covered: energies and virial (MD-Bench kernels compute forces only), the
  `-halflj` optimization, LJ-PME and switch modifiers, 2xNN/2xN/GPU kernels, MPI.

## Results

AMD Ryzen 5 7535HS (Zen 3+), GCC 12, AVX2, single precision, one thread pinned to a core,
GROMACS 2027-dev (4xM, AVX2_256), 500 iterations. Force agreement (MD-Bench vs. GROMACS):
RMS error 1.7e-6 and maximum 1.6e-5 (3.3e-5 for 24000 atoms) relative to the RMS force
(float precision; the error of either against the float64 reference is 6.6e-6 for SP).

| System | Coulomb | Kernel | us/call | computed pairs/us | useful pairs/us |
|---|---|---|---|---|---|
| 3000 atoms | Ewald | GROMACS 4xM | 1113 | 1231 | 566 |
| | | MD-Bench 4xN, half list | 1946 | 1045 | 324 |
| | | MD-Bench 4xN, full list | 2496 | 1152 | 252 |
| | RF | GROMACS 4xM | 760 | 1801 | 828 |
| | | MD-Bench 4xN, half list | 1330 | 1530 | 474 |
| 24000 atoms | Ewald | GROMACS 4xM | 9045 | 1185 | 557 |
| | | MD-Bench 4xN, half list | 13173 | 1006 | 382 |
| | RF | GROMACS 4xM | 6402 | 1675 | 787 |
| | | MD-Bench 4xN, half list | 9129 | 1452 | 552 |

Per computed pair MD-Bench's kernel runs at 85% of GROMACS' speed. The remaining
difference in time per call comes from its pair list having 51% more atom pairs.
