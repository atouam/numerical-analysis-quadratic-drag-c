#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diskfct.h"
#include "numerical_methods.h"

double m = 80.0;
double g = 9.81;
double alpha = 30.0;

double xmin = 0.0;
double xmax = 25.0;
double V0 = 0.0;
double x0 = 4000;


double acceleration(double t, double v)
{
    double tau = m / alpha;
    return g - (1.0 / tau) * v * v;
}

int velocity_analytic(const double *t, double *v, int N, double dt, double v0)
{
    double tau = m / alpha;
    double k   = sqrt(tau * g);
    double C   = (v0 - k) / (v0 + k);
    v[0] = v0;
    for (int i = 1; i < N+1; i++)
    {
        double e = exp((-2.0 * k * t[i]) / tau);
        v[i] = k * (C * e + 1.0) / (1.0 - C * e);
    }
    return 0;
}

int position_analytic(const double *t, double *x, int N, double dt)
{

    double tau = m / alpha;
    double k   = sqrt(tau * g);
    double C   = (V0 - k) / (V0 + k);

    x[0] = x0;
    for (int i = 1; i < N+1; i++)
    {
        double e = exp((-2.0 * k * t[i]) / tau);
        x[i] = x0 - k * t[i] - tau * log( (1.0 - C * e) / (1.0 - C) );
    }
    return 0;
}

int acceleration_analytic(const double *t, double *a, int N, double dt, double V0)
{
    double *v = malloc((N + 1) * sizeof *v);
    velocity_analytic(t, v, N, dt, V0);
    for (int i = 0; i < N+1; i++)
    {
        a[i] = acceleration(t[i], v[i]);
    }
    free(v);
    return 0;
}

int acceleration_forward(double *a,const double *v, int N, double dt )
{
    for (int i = 0; i < N; i++)
        a[i] = (v[i + 1] - v[i]) / dt;
    a[N] = a[N - 1];
    return 0;
}


int main(void)
{
    // Velocity
    double N_list[] = {50, 100, 200, 400};
    int nbN = (int)(sizeof N_list / sizeof N_list[0]);

    double *Errors_euler = malloc(nbN * sizeof *Errors_euler);
    double *Errors_rk2   = malloc(nbN * sizeof *Errors_rk2);

    char filename[64];

for (int i = 0; i < nbN; i++)
{
    int N = (int)N_list[i];

    double dt = 0.0;
    double *t = linspace(xmin, xmax, N, &dt);

    // --- Velocity arrays
    double *v_euler    = malloc((N + 1) * sizeof *v_euler);
    double *v_rk2      = malloc((N + 1) * sizeof *v_rk2);
    double *v_analytic = malloc((N + 1) * sizeof *v_analytic);

    // --- Position arrays
    double *x_euler    = malloc((N + 1) * sizeof *x_euler);
    double *x_rk2      = malloc((N + 1) * sizeof *x_rk2);
    double *x_analytic = malloc((N + 1) * sizeof *x_analytic);

    // --- Acceleration arrays
    double *a_euler    = malloc((N + 1) * sizeof *a_euler);      // forward diff gives N values
    double *a_rk2      = malloc((N + 1) * sizeof *a_rk2);
    double *a_analytic = malloc((N + 1) * sizeof *a_analytic);


    // Velocity
    Euler(t, v_euler, N, dt, V0, acceleration);
    RK2_Heun(t, v_rk2, N, dt, V0, acceleration);
    velocity_analytic(t, v_analytic, N, dt, V0);

    // Velocity Errors
    Errors_euler[i] = error_max_abs(v_euler, v_analytic, N);
    Errors_rk2[i]   = error_max_abs(v_rk2,   v_analytic, N);

    // Position
    Trapezoidal(x_euler, N, dt, x0, v_euler);
    Trapezoidal(x_rk2,   N, dt, x0, v_rk2);
    position_analytic(t, x_analytic, N, dt);

    // Acceleration
    acceleration_forward(a_euler, v_euler, N, dt);
    acceleration_forward(a_rk2,   v_rk2,   N, dt);
    acceleration_analytic(t, a_analytic, N, dt, V0);
    

    // Write analytical
    if (N == N_list[nbN - 1])
    {
        Write(t, a_analytic, "data/accel_analytic.dat", N);
        Write(t, x_analytic, "data/pos_analytic.dat", N);
        Write(t, v_analytic, "data/vel_analytic.dat", N);
    }

    // Write velocity
    sprintf(filename, "data/vel_euler_%d.dat", N);
    Write(t, v_euler, filename, N);

    sprintf(filename, "data/vel_rk2_%d.dat", N);
    Write(t, v_rk2, filename, N);

    // Write position
    sprintf(filename, "data/pos_euler_%d.dat", N);
    Write(t, x_euler, filename, N);

    sprintf(filename, "data/pos_rk2_%d.dat", N);
    Write(t, x_rk2, filename, N);

    // Write Acceleration
    sprintf(filename, "data/accel_euler_%d.dat", N);
    Write(t, a_euler, filename, N);

    sprintf(filename, "data/accel_rk2_%d.dat", N);
    Write(t, a_rk2, filename, N );
    

    // free
    free(t);

    free(v_euler);
    free(v_rk2);
    free(v_analytic);

    free(x_euler);
    free(x_rk2);
    free(x_analytic);

    free(a_euler);
    free(a_rk2);
    free(a_analytic);
}

    Write(N_list, Errors_euler, "data/Error_euler.dat", nbN);
    Write(N_list, Errors_rk2,   "data/Error_rk2.dat",   nbN);

    free(Errors_euler);
    free(Errors_rk2);

    // Alpha variation
    {
        int N = 400;

        double alpha_list[] = {10, 20, 30, 40};
        int nbalpha = (int)(sizeof alpha_list / sizeof alpha_list[0]);

        double *E = malloc(nbalpha * sizeof *E);

        double dt = 0.0;
        double *t = linspace(xmin, xmax, N, &dt);

        double *v = malloc((N + 1) * sizeof *v);
        double *v_analytic = malloc((N + 1) * sizeof *v_analytic);

        for (int i = 0; i < nbalpha; i++)
        {
            alpha = alpha_list[i];

            Euler(t, v, N, dt, V0, acceleration);
            velocity_analytic(t, v_analytic, N, dt, V0);

            E[i] = error_max_abs(v, v_analytic, N);
        }

        Write(alpha_list, E, "data/Error_alpha.dat", nbalpha);

        free(E);
        free(t);
        free(v);
        free(v_analytic);
    }
    
    return 0;
}