#include <iostream>
#include <cmath>
#include "kepler.cpp"
#include "numerical_methods.cpp"
#include "cauchy_problem.cpp"

using namespace std;


int main() {

    double const PI = acos(-1e0);
	int const M = 100;
	int const N = kepler_example::N;
	double period = 2 * PI;
	double tf = period * 1.5e0;

	double time_domain[M];
	for (int i = 0; i < M; i++) time_domain[i] = tf * (1e0 * i) / (1e0 * M);
    

    double u[M][N];
    for(int i = 0; i < N; i++) u[0][i] = kepler_example::u0[i];

	cauchy_problem<double>(time_domain, kepler<double>, runge_kutta<double>, u);

    cout << "Hello, World!" << endl;
    cout << "xf = (" << u[M-1][0] << ", " << u[M-1][1] << ")" << endl;
    cout << "radius = " << sqrt(pow(u[M-1][0], 2e0) + pow(u[M-1][1], 2e0)) << endl;
	getchar();

    return 0;
}