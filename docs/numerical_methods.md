# Numerical Methods

Includes several Fortran-written temporal schemes as well as
a classic Cauchy Problem interface for ODE solving.

## Cauchy

Includes the Cauchy Problem interface for ODE solving with
an standarized inputs

- (subroutine) cauchy_problem
  - (double) time_domain
  - (double function) differential_operator
  - (subroutine) temporal_scheme
  - (double) solution