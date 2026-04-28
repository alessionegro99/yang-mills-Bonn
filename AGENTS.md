# Repository Guidelines

## Project Structure & Module Organization

This is a C/Autotools codebase for lattice Yang-Mills simulations. Public
headers live in `include/`, reusable implementation code in `lib/`, and
program entry points in `src/`. Core gauge-configuration logic is split across
`lib/gauge_conf_*.c` and `include/gauge_conf.h`; group-specific routines are in
`su2*` and `sun*` files. Analysis helpers are Python scripts under `tools/`.
Autotools inputs are `configure.ac` and `Makefile.am`; generated files such as
`configure`, `Makefile.in`, and `config.h.in` are kept in the tree.

## Build, Test, and Development Commands

- `./configure`: configure a default SU(2), 4D build.
- `./configure N_c=3 ST_dim=3`: configure a different gauge group or dimension.
- `./configure --enable-use-openmp Num_threads=4`: enable OpenMP support.
- `make`: build all simulation and debug executables declared in `Makefile.am`.
- `make clean`: remove binaries and generated runtime data listed in
  `CLEANFILES`.
- `make distclean`: remove configuration output before reconfiguring from
  scratch.

If you change `configure.ac` or `Makefile.am`, regenerate the Autotools outputs
with the local Autotools toolchain before committing.

## Coding Style & Naming Conventions

The project is C99 and builds with `-Wall -Wextra -Werror -pedantic
-Wconversion`. Match nearby style: two-space indentation is common in newer
files, while older library files use brace placement with opening braces on the
next indented line. Keep names descriptive and consistent with existing
patterns, for example `Gauge_Conf`, `GParam`, `calcstaples_*`, and
`yang_mills_*`. Prefer existing macros from `include/macro.h` and compile-time
configuration from `config.h` over ad hoc constants.

## Testing Guidelines

There is no dedicated automated test suite in this checkout. Treat a clean
`make` as the baseline check, then run the most relevant debug or simulation
binary manually. Useful smoke targets include `./debug_su2`, `./debug_sun`,
`./debug_rng`, and `./conf_check`. Most simulation executables can emit a
template input file when run without an input file; use that template for small,
reproducible checks.

## Commit & Pull Request Guidelines

Recent history uses short one-line summaries, but avoid vague messages such as
`updates` or `Fix`. Prefer specific imperative summaries like `fix SU(N)
staple normalization`. Pull requests should describe the physics or numerical
behavior changed, list the configure command used, include the validation
commands and outputs checked, and mention any generated Autotools files touched.
