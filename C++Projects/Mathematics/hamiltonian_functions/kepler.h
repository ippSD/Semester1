//
// Created by imanol on 1/01/19.
//

#ifndef BODIES_KEPLER_H
#define BODIES_KEPLER_H

template  <class T>
void bind_state(T state[], T* r[2], T* v[2]);

template  <class T>
void bind_derivative(T derivative[], T* dr_dt[2], T* dv_dt[2]);

template <class T>
void kepler(T t, T y[], T yp[]);

#endif //BODIES_KEPLER_H

