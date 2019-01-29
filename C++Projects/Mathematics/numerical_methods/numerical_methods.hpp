#ifndef NUMERICAL_METHODS_HPP
#define NUMERICAL_METHODS_HPP

	template <class T>
	void runge_kutta(void (f)(T, T*, T*), T t1, T t2, T* u1, T* u2);

	template <class T>
	void euler_explicit(void (f)(T, T*, T*), T t1, T t2, T* u1, T* u2);
#endif