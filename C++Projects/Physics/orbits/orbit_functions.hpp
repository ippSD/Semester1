//
// Created by imanol on 1/01/19.
//

#ifndef ORBIT_FUNCTIONS_H
#define ORBIT_FUNCTIONS_H
#define pointer_args int bodies, T (&state)[N], int l, T (&p)[3]
#define ode_args T t, T u[N], T up[N]

// KEPLER NAMESPACE

namespace kepler_example {
	int const N = 4;
	double u0[N] = { 1e0, 0e0, 0e0, 1e0 };
}

// KEPLER POINTERS

template  <class T>
void bind_state(T state[], T* r[2], T* v[2]);

template  <class T>
void bind_derivative(T derivative[], T* dr_dt[2], T* dv_dt[2]);

// N-BODIES POINTERS

template  <class T, size_t N>
void p2m_i(T* state[N], int l, T* p);

template  <class T, size_t N>
void p2r_i(pointer_args);

template  <class T, size_t N>
void p2v_i(pointer_args);

// FUNCTIONS

template <class T, size_t N>
void kepler(ode_args);

template <class T, size_t N>
void n_bodies(ode_args);

#endif //ORBIT_FUNCTIONS_H

