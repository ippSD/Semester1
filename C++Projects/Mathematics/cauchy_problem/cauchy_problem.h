#ifndef CAUCHY_PROBLEM_H
	#define CAUCHY_PROBLEM_H
	#define F_R_RN2RN(func, tp) void(func)(tp, tp*, tp*)
	#define T_SCHEME(tp) void(temporal_scheme)(F_R_RN2RN(f,tp),tp,tp,tp*,tp*)

	template <class T, size_t M, size_t N>
	void cauchy_problem(T(&time_domain)[M], F_R_RN2RN(differential_operator, T), T_SCHEME(T), T(&solution)[M][N]);

#endif