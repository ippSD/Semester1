#include "stdafx.h"
#include <cmath>
#include "cauchy_problem.h"
#define F_R_RN2RN(func, tp) void(func)(tp, tp*, tp*)
//#define T_SCHEME(tp) void(temporal_scheme)(F_R_RN2RN(f,tp),tp,tp,tp*,tp*)
#define T_SCHEME(tp) void(temporal_scheme)(F_R_RN2RN(f,tp),tp,tp,tp*,tp*)

// void runge_kutta(void (f)(T, T*, T*), T t1, T t2, T u1[], T u2[])

template <class T, size_t M, size_t N>
void cauchy_problem(T (&time_domain)[M], F_R_RN2RN(differential_operator, T), T_SCHEME(T), T (&solution)[M][N]) {

	for (int i = 0; i < M - 1; i++) {
		temporal_scheme(differential_operator, time_domain[i], time_domain[i + 1], solution[i], solution[i + 1]);
	}
}