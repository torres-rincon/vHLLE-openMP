# Running SMASH-vHLLE with OpenMP-enabled vHLLE

This note explains how to run the SMASH-vHLLE hybrid setup when the vHLLE
hydrodynamic evolution has been compiled with OpenMP support.

The SMASH commands and SMASH configuration files are unchanged. The only new
control is the number of OpenMP threads used by `hlle_visc`.

## 1. Compile vHLLE with OpenMP

From the vHLLE repository:

```sh
make OPENMP=1
```

`OPENMP=1` only enables OpenMP support at compile time. It does not select the
number of threads used at runtime.

To rebuild from scratch:

```sh
make clean
make OPENMP=1
```

On Linux this uses the compiler OpenMP flags directly. On macOS, the Makefile
expects an OpenMP runtime such as Homebrew `libomp` when using Apple Clang.

## 2. Configure the hybrid CMake build

Create and configure a hybrid build directory as usual, but pass
`VHLLE_OMP_NUM_THREADS` to choose how many threads vHLLE should use during the
hydro step:

```sh
cmake -S . -B build-hybrid \
  -DSMASH_PATH=/path/to/smash/build \
  -DVHLLE_PATH=/path/to/vhlle \
  -DVHLLE_PARAMS_PATH=/path/to/vhlle \
  -DSAMPLER_PATH=/path/to/sampler/build \
  -DVHLLE_OMP_NUM_THREADS=8
```

With this configuration, every CMake-launched `hlle_visc` hydro command runs
with:

```sh
OMP_NUM_THREADS=8
OMP_DYNAMIC=FALSE
```

`OMP_DYNAMIC=FALSE` prevents the OpenMP runtime from dynamically changing the
number of threads, so `VHLLE_OMP_NUM_THREADS=8` really means that vHLLE should
use 8 threads unless the system imposes a stricter limit.

If `VHLLE_OMP_NUM_THREADS` is omitted, CMake launches `hlle_visc` as before and
the OpenMP runtime uses its default behavior or the external environment.

## 3. Run the custom hybrid targets

For the custom SMASH initial-condition workflow, the relevant targets are:

```sh
cmake --build build-hybrid --target custom_IC
cmake --build build-hybrid --target custom_hydro
cmake --build build-hybrid --target custom_sampler
cmake --build build-hybrid --target custom_afterburner
```

To run the hydro stage with 8 OpenMP threads per vHLLE process:

```sh
cmake --build build-hybrid --target custom_hydro --parallel 1
```

The `--parallel` option controls how many build jobs CMake runs at the same
time. It is separate from `VHLLE_OMP_NUM_THREADS`. For example, running
`--parallel 4` with `VHLLE_OMP_NUM_THREADS=8` may launch up to four hydro jobs
at once, each using eight OpenMP threads.

For long production runs, choose these values so that:

```text
number of simultaneous CMake jobs × VHLLE_OMP_NUM_THREADS <= available CPU cores
```

For example, on an 8-core node:

```sh
cmake --build build-hybrid --target custom_hydro --parallel 1
```

with:

```sh
-DVHLLE_OMP_NUM_THREADS=8
```

uses one hydro event at a time with eight vHLLE threads.

## 4. Changing the number of vHLLE threads

To switch from 8 to another number of threads, reconfigure the same build
directory:

```sh
cmake -S . -B build-hybrid -DVHLLE_OMP_NUM_THREADS=16
```

Then run the desired target again:

```sh
cmake --build build-hybrid --target custom_hydro --parallel 1
```

No SMASH configuration needs to be changed.
