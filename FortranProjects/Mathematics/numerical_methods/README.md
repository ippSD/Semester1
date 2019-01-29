# numerical_methods

Includes numerical methods for solving ODE with the following structure:

<p style="text-align: center;"><b> ∫du = ∫F(u, t)·dt</b></p>
<p style="text-align: center;"><b>u(0) = u0</b></p>

* Cauchy Problem: Solves the initial value problem shown above.
* Richardson Extrapolation: solves the Cauchy Problem for a more refined
time domain in order to estimate the integration error.