#include <stdlib.h>
#include <math.h>

double *linspace(double xmin, double xmax, int N, double *dt)
{
    double h = (xmax - xmin) / N;
    if (dt) *dt = h;

    double *X = malloc((N + 1) * sizeof *X);
    if (!X) return NULL;

    for (int i = 0; i <= N; i++)
        X[i] = xmin + i * h;

    return X;
}


int Euler(const double *t, double *Y, int N, double dt, double Y0, double (*f)(double, double))
{
    Y[0] = Y0;
    for (int i = 0; i < N; i++)
        Y[i + 1] = Y[i] + dt * f(t[i], Y[i]);
    return 0;
}

int RK2_Heun(const double *t, double *Y, int N, double dt, double Y0, double (*f)(double, double))
{
    Y[0] = Y0;
    for (int i = 0; i < N; i++)
    {
        double ti = t[i];
        double yi = Y[i];

        double k1 = f(ti, yi);
        double k2 = f(ti + 0.5 * dt, yi + 0.5 * dt * k1);

        Y[i + 1] = yi + dt * k2;
    }
    return 0;
}


int Trapezoidal( double *x, int N, double dt, double x0, const double *v)
{
    x[0] = x0;
    for (int i = 0; i < N; i++)
        x[i + 1] = x[i] - 0.5 * dt * (v[i] + v[i + 1]);
    return 0;
}


double error_max_abs(const double *Y_num, const double *Y_ref, int N)
{
    double emax = 0.0;
    for (int i = 0; i <= N; i++)
    {
        double e = fabs(Y_num[i] - Y_ref[i]);
        if (e > emax) emax = e;
    }
    return emax;
}