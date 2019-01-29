#pragma once

namespace NumericalMethods
{
	template <class T>
	void runge_kutta(void (f)(T, T*, T*), T t1, T t2, T* u1, T* u2);

	template <class T>
	void euler_explicit(void (f)(T, T*, T*), T t1, T t2, T* u1, T* u2);
}