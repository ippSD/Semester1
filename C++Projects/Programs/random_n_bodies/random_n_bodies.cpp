// random_n_bodies.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "pch.h"
#include <iostream>
#include <cmath>
#include <time.h>
#include "orbit_functions.cpp"
#include "numerical_methods.cpp"
#include "cauchy_problem.cpp"

#define loop(i, n) for(int i = 0; i < n; i++)

using namespace std;

int main()
{
	// Total Number of Bodies.
	int const BODIES = 800;

	// Total time steps + Initial time.
	int const M = 101;

	// State vector size.
	const size_t N = BODIES * 7;

	// Final time.
	double const TF = 1e3;

	// Time Domain.
	double time_domain[M];

	// InitialConditions as vector.
	double u0[N];

	// U0 pointer;
	double* pu0[N];

	// Position pointer.
	double* r[3];

	// Velocity pointer.
	double* v[3];

	double* x[3];

	// State Matrix.
	double u[M][N];

	loop(i, N) pu0[i] = &u0[i];
	loop(i, M) time_domain[i] = i * TF / (M - 1e0);


	//Random Seed based on time;
	srand(time(0));
	//Get random initial values in[0, 1)
	loop(i, N) u0[i] = rand();
	//Get negative values
	loop(i, N) u0[i] = (u0[i] - 5e-1) * 2e0;
	//Positive mass with a classic factor of 1E6 km ^ 3 s^-2
	loop(i, BODIES) u0[i] = abs(u0[i]) * 1e6;
	//Positions of order 1E4 km
	loop(i, 3 * BODIES) u0[BODIES + i] *= 1e4;
	//Velocities of order 1E0 km / s
	loop(i, 3 * BODIES) u0[4 * BODIES + i] *= 5e-1;

	//Z coordinates are null (2D)
	loop(i, BODIES) {
		p2r_i(BODIES, pu0, i, r);
		*r[2] = 0e0;
		p2v_i(BODIES, pu0, i, v);
		*v[2] = 0e0;
	}

	//Solve
	loop(i, N) cout << u0[i] << endl;
	cauchy_problem<double>(time_domain, n_bodies<double, N>, runge_kutta<double>, u);
}