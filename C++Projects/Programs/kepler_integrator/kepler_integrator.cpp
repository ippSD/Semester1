#include <iostream>
#include <cmath>
//#include "orbit_functions.hpp"
//#include "numerical_methods.h"
//#include "cauchy_problem.h"
#include "cauchy_problem.cpp"
#include "orbit_functions.cpp"
#include "numerical_methods.cpp"

using namespace std;

int main() {

    double const PI = acos(-1e0);
	int const M = 100;
	int const N = kepler_example::N;
	double period = 2 * PI;
	double tf = period * 1.5e0;
	double time_domain[M];
	double u_rk[M][N], u_ee[M][N];
	double* x[M];
	double* y[M];

	for (int i = 0; i < M; i++) time_domain[i] = tf * (1e0 * i) / (1e0 * M);

	getchar();
   
	for (int i = 0; i < N; i++) {
		u_rk[0][i] = kepler_example::u0[i];
		u_ee[0][i] = kepler_example::u0[i];
	}
	for (int i = 0; i < M; i++) {
		x[i] = &u_rk[i][0];
		y[i] = &u_rk[i][1];
	}

	getchar();

	cauchy_problem<double>(time_domain, kepler<double, N>, runge_kutta<double>, u_rk);

	getchar();

	cauchy_problem<double>(time_domain, kepler<double, N>, euler_explicit<double>, u_ee);

    cout << "Kepler Integrator" << endl;
	cout << "Position and radius after 1 orbit period." << endl << endl;

	cout << "-> Runge-Kutta 4" << endl;
    cout << "xf = (" << u_rk[M-1][0] << ", " << u_rk[M-1][1] << ")" << endl;
    cout << "radius = " << sqrt(pow(u_rk[M-1][0], 2e0) + pow(u_rk[M-1][1], 2e0)) << endl;

	cout << "-> Explicit Euler" << endl;
	cout << "xf = (" << u_ee[M - 1][0] << ", " << u_ee[M - 1][1] << ")" << endl;
	cout << "radius = " << sqrt(pow(u_ee[M - 1][0], 2e0) + pow(u_ee[M - 1][1], 2e0)) << endl;


	getchar();

    return 0;
}