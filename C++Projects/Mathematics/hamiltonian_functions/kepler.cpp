//
// Created by imanol on 1/01/19.
//
#include "stdafx.h"
#include "kepler.h"

namespace kepler_example {
	int const N = 4;
	double u0[N] = { 1e0, 0e0, 0e0, 1e0 };
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

template <class T>
void kepler(T t, T y[], T yp[]) {
	T* r[2];
	T* v[2];
	T* dr_dt[2];
	T* dv_dt[2];
	bind_state(y, r, v);
	bind_derivative(yp, dr_dt, dv_dt);

	T d2 = pow(*r[0], 2e0) + pow(*r[1], 2e0);
	*dr_dt[0] = *v[0];
	*dr_dt[1] = *v[1];
	*dv_dt[0] = -*r[0] / d2;
	*dv_dt[1] = -*r[1] / d2;
};
