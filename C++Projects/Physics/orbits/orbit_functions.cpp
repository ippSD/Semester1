//
// Created by imanol on 1/01/19.
//
#include "stdafx.h"
#include "orbit_functions.hpp"
#include <math.h>
#include <iostream>

using namespace std;

template  <class T>
T norm2(T r[3]) {
	return sqrt(pow(r[0], 2e0) + pow(r[1], 2e0) + pow(r[2], 2e0));
}

template  <class T>
T norm2(T r1[3], T r2[3]) {
	return sqrt(pow(r1[0] - r2[0], 2e0) + pow(r1[1] - r2[1], 2e0) + pow(r1[2] - r2[2], 2e0));
}

template  <class T>
void bind_state(T state[], T* r[2], T* v[2]) {
	r[0] = &state[0];
	r[1] = &state[1];
	v[0] = &state[2];
	v[1] = &state[3];
}

template  <class T>
void bind_derivative(T derivative[], T* dr_dt[2], T* dv_dt[2]) {
	dr_dt[0] = &derivative[0];
	dr_dt[1] = &derivative[1];
	dv_dt[0] = &derivative[2];
	dv_dt[1] = &derivative[3];
}

template <class T, size_t N>
void kepler(ode_args) {
	T* r[2];
	T* v[2];
	T* dr_dt[2];
	T* dv_dt[2];
	bind_state(u, r, v);
	bind_derivative(up, dr_dt, dv_dt);

	T d2 = pow(*r[0], 2e0) + pow(*r[1], 2e0);
	*dr_dt[0] = *v[0];
	*dr_dt[1] = *v[1];
	*dv_dt[0] = -*r[0] / d2;
	*dv_dt[1] = -*r[1] / d2;
};

/*
N BODIES PROBLEM
*/

template  <class T, size_t N>
void p2m_i(T (&state)[N], int l, T (&p)) {
	p = state[l];
}

template  <class T, size_t N>
void p2r_i(pointer_args) {
	for (int i = 0; i < 3; i++) p[i] = state[bodies + 3*l + i];
}

template  <class T, size_t N>
void p2v_i(pointer_args) {
	for (int i = 0; i < 3; i++) p[i] = state[4*bodies + 3*l + i];
}

template <class T,  size_t N>
void n_bodies(ode_args) {
	T *m_j;
	T *r_i[3], *r_j[3], *v_i[3];
	T *dr_dt_i[3], *dv_dt_i[3];
	T d, a_j;
	T *pu[N];
	T *pup[N];
	for (int i = 0; i < N; i++) {
		pu[i] = &u[i];
		pup[i] = &up[i];
	}
	int bodies = N / 7;
	
	for (int i = 0; i < bodies; i++)
	{
		p2r_i(bodies, pu,  i,     r_i);  // Pointer to i-BODY POSITION
		p2r_i(bodies, pup, i, dr_dt_i);  // Pointer to i-BODY POSITION DERIVATIVE
		p2v_i(bodies, pu,  i,     v_i);  // Pointer to i-BODY VELOCITY
		p2v_i(bodies, pup, i, dv_dt_i);  // Pointer to i-BODY VELOCITY DERIVATIVE

		*dr_dt_i = *v_i;
		for (int k = 0; k < 3; k++) *dv_dt_i[k] = 0e0;

		for (int j = 0; j < bodies; j++)
		{
			p2m_i(pu, j, m_j);  // Pointer to j-BODY MASS
			p2r_i(bodies, pu, j, r_j);  // Pointer to j-BODY POSITION
			if (j != i)
			{
				d = norm2(*r_i, *r_j);
				for (int k = 0; k < 3; k++)
				{
					a_j = *m_j;
					a_j *= (*r_i[k] - *r_j[k]);
					a_j /= pow(d, 2e0);
					*dv_dt_i[k] -= a_j;
				}
				
			}
		}
	}
};