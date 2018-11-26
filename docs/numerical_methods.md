# Numerical Methods

Includes several Fortran-written temporal schemes as well as
a classic Cauchy Problem interface for ODE solving.

## Cauchy

Includes the Cauchy Problem interface for ODE solving with
an standarized inputs

- (subroutine) cauchy_problem(time_domain, differential_operator, temporal_scheme, solution)
  -  (double precision, intent(in)) time_domain(0:) : 'M'-sized sorted time partition. When the current iteration time approaches the current state vector is saved.
  -  (procedure(odes)) differential_operator : derivative function of the state vector which must match 'ode' function template.
  -  (procedure(scheme)) temporal_scheme : numerical scheme to be used for temporal integration which must match 'scheme' subroutine template.
  -  (double precision, intent(inout)) solution(0:,:) : 'M'x'N'-sized state matrix. The 'N'-sized state vector's 'm'-th step is stored in it. Index 0 corresponds to the input initial condition.

## Differential Operator

Derivative function of the state vector matching the following template:

- (function) odes(u, t) result(f)
  - (double precision, intent(in)) u(:) : Current iteration state vector.
  - (double precision, intent(in)) t : Current iteration time.
  - (double precision) f(:) : Derivate of the state vector with same size.

## Temporal Scheme

Temporal integration numerical scheme which matches the following template:

- (subroutine) scheme(f, t1, t2, u1, u2)
  - (procedure(odes)) f : Same as 'Differential Operator'.
  - (double precision, intent(in)) t1 : Current iteration time.
  - (double precision, intent(in)) t2 : Time when the solution will be given.
  - (double precision, intent(in)) u1(:) : Current state vector.
  - (double precision, intent(out)) u2(:) : Solution at time 't2'. Same size as 'u1'.
