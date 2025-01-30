#ifndef MYNUMERICS_H
#define MYNUMERICS_H
#include <gsl/gsl_rng.h>

typedef int ODE_FUNC(double,const double[],double[],void*);
typedef int SDEfunc(double, const double[], double[], double[], void*);
typedef struct{
    double x;
    double y;
}Tuple;
typedef struct{
    int size; //anzahl von Bins
    double low; //unterste Intervallgrenze
    double high; //oberste Intervallgrenze
    double* bin; //Array mit Klassenhäufigkeiten
} Hist;


//Gaussian integrand
double mn_gaussian(double x);

//Cosinus hyperbolicus
double mn_cosh(double x);

//Signum function
double mn_sign(double x);

//Euklidische Norm von vecotr (array Form)
double mn_euclid_norm(int dim, double *x);

//Factorial natürlicher Zahl
int mn_int_factorial(int n);

//Integrand von Gammafunktion, x reell
double mn_gamma_integrand(double t, double x);

//Midpoint Integrator 1D
double mn_integrate_midpoint(double left, double right, double N, double integrand(double y));

//Simpsonregel Integrator 1D
double mn_integrate_simpson(double left, double right, double N, double integrand(double y));

//Simpsonregel Integrator für 2D
double mn_integrate_simpson_two(double left, double right, double delta_x, double integrand(double, double), double k);

//Gammafunktion für reele Argumente
double mn_gamma_function_real(double x, int max_iter, double rel_tol);

//Gaussian error function Midpoint und Simpsonregel
double mn_erf_midpoint(double x, double delta_x);
double mn_erf_simpson(double x, double delta_x);

//cosh Integral mit Midpoint und Simpson
//int cosh(x)dx from 0 to 1
double mn_cosh_midpoint(double delta_x);
double mn_cosh_simpson(double delta_x);

//numerical differentiator
double mn_diff(double x, double delta, double func(double));

//forward difference derivative
double mn_frwrd_diff(double x, double delta, double func(double));

//box function = 1 for x in [-0.5,0.5] else 0
double mn_rect(double x);

//Newton-Raphson
double mn_find_root(double func(double), double x0, double delta, const double rel_tol, const int max_iter);

//Euler stepper
void mn_euler_step(double t, double dt, double y[], double dy[], ODE_FUNC func, int dim, void* params);

//RK2 stepper
void mn_rk2_step(double t, double delta_t, double y[], double dy[], ODE_FUNC func, int dim, void *params);

//RK4 stepper
void mn_rk4_step(double t, double delta_t, double y[], double dy[], ODE_FUNC func, int dim, void *params);

//VV stepper
void mn_velocity_verlet_step(double t, double delta_t, double y[], double dy[], ODE_FUNC func, int dim, void *params);

//Zwei gleichverteilte Zufallszahlen in (-1,1)
Tuple mn_random_gaussian(gsl_rng* generator);

//Euler-Maruyama
void mn_euler_maruyama_step(double t, double delta_t, double y[], SDEfunc func, int dimension, gsl_rng* gen, void *params);

/*-------------------------------------------- Histogramme ---------------------------------------------------------*/

//Speicher allokieren
Hist* hist_malloc(int size, double low, double high);

//Speicher freigeben
void hist_free(Hist* H);

//Histogramm mit Null füllen
void hist_reset(Hist* H);

//Wert addieren
void hist_add(Hist* H, double x);

//Histogramm in CSV Datei schreib
void hist_fprintf(Hist* H, FILE* file);

//Summe aller Samples
double hist_sum(Hist* H);

//wie hist_fprintf nur mir Relativen Dichten, statt absoluten
void hist_fprintf_rel(Hist* H, FILE* file);

//Erwartungswert von Histogramm
double hist_mean(Hist* H);

//Standardabweichung von Histogramm
double hist_std(Hist* H);

#endif