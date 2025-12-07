# CHANGELOG

All notable changes to `dae-cpp` project will be documented in this file.

## [`develop` branch -- 2.3.0]

### Added

- Option `linear_system_scaling` (can be `true` or `false`) to rescale linear system matrix for better convergence (`false` by default).
- Linear system matrix (row) scaling. Matrix scaling can improve stability of the linear solver in some cases but comes with a slight performance penalty.
- Linear system matrix prune if scaling is enabled.

### Fixed

- The solver cannot reach relative tolerance `rtol` in some cases and stops with the "solution diverged" error (fixed by updating internal tolerances).
- Incompatibility of `autodiff` library with the latest version of Eigen (fixed by altering `autodiff`).
- Typo in `daecpp::solver_command::stop_integration` enum (`stop_intergration` -> `stop_integration`).

### Changed

- Eigen to version 5.0.0.
- Renamed `TESTING` macro definition to `DAECPP_TESTING` to avoid potential clash.
- Linear system matrix pattern now analysed only once at the first iteration.
- Pre-allocate vector of `dual` numbers in `JacobianMatrixShape` class to improve performance of the Jacobian computed from the user-defined shape.
- Updated internal tolerances used in the solver for the convergence check against relative tolerance `rtol`.

### Removed

- Conversion from Eigen to dae-cpp matrix format and back in automatic Jacobian class (an attempt to speed it up).

## [2.2.0]

Current stable version.

[CHANGELOG](https://dae-cpp.github.io/CHANGELOG.html)
