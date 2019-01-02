#include "stdafx.h"
#include "rkf45.hpp"


template <class T>
void runge_kutta(void (f)(T, T*, T*), T t1, T t2, T* u1, T* u2) {
	int const m = *(&u1 + 1) - u1;
	T t = t1;
	T* du_dt = new T[m];
	T rel_err = 1e-5;
	T const ABS_ERR = 1e-5;
	for (int i = 0; i < m; i++) u2[i] = u1[i];
	r8_rkf45(f, m, u2, du_dt, &t, t2, &rel_err, ABS_ERR, 1);
}

template <>
void runge_kutta<float>(void (f)(float, float*, float*), float t1, float t2, float* u1, float* u2) {
	int const m = *(&u1 + 1) - u1;
	float t = t1;
	float* du_dt = new float[m];
	float rel_err = 1e-5;
	float const ABS_ERR = 1e-5;
	for (int i = 0; i < m; i++) u2[i] = u1[i];
	r4_rkf45(f, m, u2, du_dt, &t, t2, &rel_err, ABS_ERR, 1);
}