template<typename T>
T r4_abs ( T x );
template<typename T>
T r4_epsilon ( );
template<typename T>
void r4_fehl ( void f ( T t, T y[], T yp[] ), int neqn,
  T y[], T t, T h, T yp[], T f1[], T f2[], T f3[], 
  T f4[], T f5[], T s[] );
template<typename T>
T r4_max ( T x, T y );
template<typename T>
T r4_min ( T x, T y );
template<typename T>
int r4_rkf45 ( void f ( T t, T y[], T yp[] ), int neqn,
  T y[], T yp[], T *t, T tout, T *relerr, T abserr, 
  int flag );
template<typename T>
T r4_sign ( T x );

double r8_abs ( double x );
double r8_epsilon ( );
void r8_fehl ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double t, double h, double yp[], double f1[], double f2[], double f3[], 
  double f4[], double f5[], double s[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_rkf45 ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double yp[], double *t, double tout, double *relerr, double abserr, 
  int flag );
double r8_sign ( double x );

void timestamp ( );
