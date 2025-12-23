#ifndef NUMERICAL_H
#define NUMERICAL_H

// Time grid
double *linspace(double xmin, double xmax, int N, double *dt);

// ODE solvers
int Euler(const double *t, double *Y, int N, double dt, double Y0, double (*f)(double, double));
int RK2_Heun(const double *t, double *Y, int N, double dt, double Y0, double (*f)(double, double));

// Integration methods
int Trapezoidal(double *x, int N, double dt, double x0, const double *v);

// Error calculation
double error_max_abs(const double *Y_num, const double *Y_ref, int N);
#endif 