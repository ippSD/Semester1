#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H
	template <class T>
	void runge_kutta(void (f)(T, T*, T*), T t1, T t2, T* u1, T* u2);
#endif