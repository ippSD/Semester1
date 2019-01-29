# temporal_schemes

Project including several ODE propagators with the following form:

<p style="text-align: center;"><b>subroutine temporal_scheme(f, t1, t2, u1, u2)</b></p>

* <b>f</b>: Derivative function.
* **t1**: Time of the current state vector, **u1**.
* **t2**: Evaluation time of the state vector at the next step, **u2**.
* **u1**: State vector at **t1**.
* **u2**: State vector at the next step.

Included propagators:
* Explicit Euler
* Inverse Euler
* Runge Kutta 4