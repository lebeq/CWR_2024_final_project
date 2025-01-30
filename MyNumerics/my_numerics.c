#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "my_numerics.h"

typedef int ODE_FUNC(double,const double[],double[],void*);
typedef int SDEfunc(double, const double[], double[], double[], void*);


//Integrand für ERF
double mn_gaussian(double x){
    return exp(-x*x);
}
//Integrand für cosh
double mn_cosh(double x){
    return 0.5*(exp(x) + exp(-x));
}

//Factorial natürlicher Zahl
int mn_int_factorial(int n){
    int result = 1;
    if(n == 0){
        return result;
    }
    for(int i = n; i>0; i--){
        result *= i;
    }
    printf("factorial of %d is %d\n", n, result);
    return result;
}

//Euklidische Norm von einem vector (array)
double mn_euclid_norm(int dim, double *x){
    double res = 0.0;
    for(int i = 0; i<dim; i++){
        res += x[i]*x[i];
    }
    return sqrt(res);
}

//signum function
double mn_sign(double x){
    if(x == 0.0){
        return 0.0;
    }
    if(x < 0.0){
        return -1.0;
    }
    if(x > 0.0){
        return 1.0;
    }
    return EXIT_FAILURE;
}

//Integration mit midpoint regel
double mn_integrate_midpoint(double left, double right, double N, double integrand(double y)){
    
    //check ob Grenzen richtig sind
    if(left > right){
        return -1*mn_integrate_midpoint(right,left,N,integrand);
    }
    double delta_x = (left-right)/N;
    double sum = 0.0;

    for(int k = 0; k < N; k++){
        double y = left + (k-0.5)*delta_x;
        double f = integrand(y);
        double A = f*delta_x;
        sum += A;
    }
    return sum;
}

//integration mit simpsonregel
double mn_integrate_simpson(double left, double right, double N, double integrand(double y)){
    
    //check ob Grenzen richtig sind
    if(left > right){
        return -1*mn_integrate_simpson(right,left,N,integrand);
    }

    double delta_x = fabs(right-left)/N;
    double sum = 0.0;

    for(int k = 0; k<N; k++){
        double x1 = left + k*delta_x;
        double m = x1 + 0.5*delta_x;
        double x2 = x1 + delta_x;
        double f = integrand(x1) + 4*integrand(m) + integrand(x2);
        double A = f*delta_x/6.0;
        sum += A;
    }
    return sum;
}
//integration mi simpsonregel für Integranden mit 2 Args
//Für Aufgabe 7
double mn_integrate_simpson_two(double left, double right, double delta_x, double integrand(double, double), double k){
    double sum = 0.0;
    //check ob Grenzen richitg sind
    if(left > right){
        return -1*mn_integrate_simpson_two(right,left,delta_x,integrand,k);
    }
    int N = (right-left)/delta_x;

    for(int i = 0; i<N; i++){
        double x1 = left + i*delta_x;
        double m = x1 + 0.5*delta_x;
        double x2 = x1 + delta_x;
        double f = integrand(x1,k) + 4.0*integrand(m,k) + integrand(x2,k);
        double A = f*delta_x/6.0;
        sum += A;
    }
    return sum;
}

//Integrand für Gamm Fkt reel
double mn_gamma_integrand(double t, double x){
    //t is the variable of integration, x is argument of Gamma
    return pow(t,x-1)*exp(-t);
}

//Gmmafunktion für reele Arg
double mn_gamma_function_real(double x, int max_iter, double rel_tol){
    double right = 1e2;
    double old = 0.0;
    double res = mn_integrate_simpson_two(0.0,right,1e-3,mn_gamma_integrand,x);
    int iter = 0;
    double delta_right = 1e3;
    while(iter < max_iter || fabs(old - res) > rel_tol){
        old = res;
        right += delta_right;
        res = mn_integrate_simpson_two(0.0,right,1e-3,mn_gamma_integrand,x);
        iter++;
    }
    return res;
}

//ERF mit midpoint
double mn_erf_midpoint(double x, double delta_x){
    double c = 2*pow(sqrt(M_PI),-1);
    double N = fabs(x)/delta_x;
    
    double ERF = mn_integrate_midpoint(0.0,x,N,mn_gaussian);

    return c*ERF;
}

//ERF mit Simpson
double mn_erf_simpson(double x, double delta_x){
    double c = 2*pow(sqrt(M_PI),-1);
    double N = fabs(x)/delta_x;

    double ERF = mn_integrate_simpson(0.0,x,N,mn_gaussian);

    return c*ERF;
}

//cosh mit midpoint
double mn_cosh_midpoint(double delta_x){
    double N = 1/delta_x;
    
    double cosh = mn_integrate_midpoint(0.0,1.0,N,mn_cosh);

    return cosh;
}

double mn_cosh_simpson(double delta_x){
    double N = 1/delta_x;

    double cosh = mn_integrate_simpson(0.0,1.0,N,mn_cosh);

    return cosh;
}

//Trägt den abs Fehler zwisch ERF und cosh integral abh von delta_x in CSV ein
//x-achse: delta_x
//y-achse: abs(erf-cosh)
//midpoint fehler
int mn_err_plot_mid(int max_del, FILE* data){
    double step = (max_del-1.0)/1e5;
    for(double i = 1.0; i < max_del; i += step){
        double erf = mn_erf_midpoint(1.0,i);
        double cosh = mn_cosh_midpoint(i);
        double abs_diff = fabs(erf-cosh);
        fprintf(data, "%.16f, %.16f\n", i, abs_diff);
    }
    return 0;
}

//simpsonregel fehler
int mn_err_plot_simp(int max_del, FILE* data){
    double step = (max_del-1.0)/1e5;
    for(double i = 1.0; i < max_del; i += step){
        double erf = mn_erf_simpson(1.0,i);
        double cosh = mn_cosh_simpson(i);
        double abs_diff = fabs(erf-cosh);
        fprintf(data, "%.16f, %.16f\n", i, abs_diff);
    }
    return 0;
}


//Numerische Ableitung
double mn_diff(double x, double delta, double func(double)){
    return (func(x+delta)-func(x-delta))/(2.0*delta);
}

//Vorwärtsdifferenz Ablt
double mn_frwrd_diff(double x, double delta, double func(double)){
    return (func(x+delta)-func(x))/delta;
}

//rect fkt: if |x|<= 0.5 rect = 1, else rect = 0
double mn_rect(double x){
    if(fabs(x)<= 0.5){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

//Newton-Raphson
double mn_find_root(double func(double), double x0, double delta, const double rel_tol, const int max_iter){
    double x = x0;
    int iteration = 0;

    while(iteration < max_iter){
        //auswertung von Fkt und Ableitung
        double f = func(x);
        double df = mn_diff(x,delta,func);

        //Speichert alten wert für toleranzvergl.
        double x_old = x;

        //Iterationsschritt
        x -= f/df;

        //checkt Toleranz
        if(fabs(x_old-x)/fabs(x) < rel_tol){
            break;
        }

        iteration++;
    }
    return x;
}

/* 
------- ANMERKUNG ZU ODE STEPPERN: ------------
Bei beiden Steppern wird dy[] übergeben, damit es einhetlich ist und die Simulationsvprgänge
ggf. als alleinestehende Methode außerhalb von main geschrieben werden können
*/

//Euler Verfahren für ODE 1. Grades
void mn_euler_step(double t, double dt, double y[], double dy[], ODE_FUNC func, int dim, void* params){
    func(t,y,dy,params);
    for(int i = 0; i < dim; i++){
        y[i] = y[i] + dy[i]*dt;
    }
}

//Runge-Kutta 2. Ordnung (RK2)
void mn_rk2_step(double t, double delta_t, double y[], double dy[], ODE_FUNC func, int dim, void *params){
    //Speicher für Stützstellen
    double *k_1 = malloc(dim*sizeof(double));
    double *dy_1 = malloc(dim*sizeof(double));

    //k_1 und output
    func(t,y,dy_1,params);
    for(int i = 0; i<dim; i++){
        k_1[i] = delta_t * dy_1[i];
    }

    //y + k_1/2 berechnen
    double *temp_y = malloc(dim*sizeof(double));
    for(int i = 0; i<dim; i++){
        temp_y[i] = y[i] + k_1[i]/2.0;    
    }

    //neues y berechnen
    //reused dy_1 als output für ODE func
    func(t + delta_t/2.0,temp_y,dy_1,params);
    for(int i = 0; i<dim; i++){
        y[i] = y[i] + delta_t*dy_1[i];
    }

    //um das Warning beim kompilieren zu sparen
    dy[0] = 1.0*dy[0];

    free(k_1);
    free(dy_1);
    free(temp_y);
}


//Runge-Kutta 4. Ordnung (RK4)
void mn_rk4_step(double t, double delta_t, double y[], double dy[], ODE_FUNC func, int dim, void *params){
    double size = dim*sizeof(double);
    //Stützstellen initialisieren
    double *k_1 = malloc(size);
    double *k_2 = malloc(size);
    double *k_3 = malloc(size);
    double *k_4 = malloc(size);
    
    //Arrays für die Ausgaben von func für die k_i's
    double *dy_1 = malloc(size);
    double *dy_2 = malloc(size);
    double *dy_3 = malloc(size);
    double *dy_4 = malloc(size);

    //Stützstellen berechenen
    //k_1
    func(t,y,dy_1,params);
    for(int i = 0; i<dim; i++){
        k_1[i] = delta_t*dy_1[i];
    }

    //k_2
    //y+k1/2, t+dt/2
    double *temp_k1 = malloc(size);
    double *temp_y1 = malloc(size);
    for(int i = 0; i < dim; i++){
        temp_k1[i] = k_1[i]/2.0;
        temp_y1[i] = y[i] + temp_k1[i];
    }
    func(t+delta_t/2.0,temp_y1,dy_2,params);
    for(int i = 0; i < dim; i++){
        k_2[i] = delta_t*dy_2[i];
    }

    //k_3
    //y+k2/2, t+dt/2
    double *temp_k3 = malloc(size);
    double *temp_y3 = malloc(size);
    for(int i = 0; i < dim; i++){
        temp_k3[i] = k_2[i]/2.0;
        temp_y3[i] = y[i] + temp_k3[i];
    }
    func(t+delta_t/2.0,temp_y3,dy_3,params);
    for(int i = 0; i < dim; i++){
        k_3[i] = delta_t*dy_3[i];
    }

    //k_4
    //y+k3, t+dt
    double *temp_y4 = malloc(size);
    for(int i = 0; i < dim; i++){
        temp_y4[i] = y[i] + k_3[i];
    }
    func(t+delta_t,temp_y4,dy_4,params);
    for(int i = 0; i < dim; i++){
        k_4[i] = delta_t*dy_4[i];
    }

    // RK4 Step
    for(int i = 0; i < dim; i++){
        double sum = k_1[i] + 2.0*k_2[i] + 2.0*k_3[i] + k_4[i];
        //STEP
        y[i] = y[i] + (1.0/6.0)*sum;
    }

    //um das Warning beim kompilieren zu sparen
    dy[0] = 1.0*dy[0];

    free(k_1);
    free(k_2);
    free(k_3);
    free(k_4);
    free(dy_1);
    free(dy_2);
    free(dy_3);
    free(dy_4);
    free(temp_k1);
    free(temp_k3);
    free(temp_y1);
    free(temp_y3);
    free(temp_y4);
}

//Velocity-Verlet Integrator (VV)
void mn_velocity_verlet_step(double t, double delta_t, double y[], double dy[], ODE_FUNC func, int dim, void *params){
    //bis N gehen die Positionen in y, ab N die Geschw.
    int N = dim/2;
    //Speicher allokieren
    double *a = malloc(dim*sizeof(double));
    double *a_cp = malloc(dim*sizeof(double));

    //positionen berechnen
    func(t,y,a,params);
    //die Beschleunigungen stehen in der zweiten Hälfte von a[]
    //berechene von x_{i+1} und a_i kopieren
    for(int i = 0; i<N; i++){
        y[i] = y[i] + y[N+i]*delta_t + a[N+i]*(delta_t*delta_t)/2.0;
        //a_i in a_cp speichern
        a_cp[i] = a[i];
        a_cp[N+i] = a[N+i];
    }
    //a_{i+1} und v_{i+1} berechenen
    func(t+delta_t,y,a,params);
    for(int i = 0; i<N; i++){
        y[N+i] = y[N+i] + (a[N+i]+a_cp[N+i])*(delta_t/2.0);
    }

    //um das Warning beim kompilieren zu sparen
    dy[0] = 1.0*dy[0];
    
    free(a);
    free(a_cp);
}

//Polarmethode, gibt zwei normalverteilte Zufallszahlen zurück
Tuple mn_random_gaussian(gsl_rng* generator){
    Tuple result;

    //get two uniform random numbers in (-1,1)
    double u = gsl_rng_uniform(generator)*2.0 - 1.0;
    double v = gsl_rng_uniform(generator)*2.0 - 1.0;
    double vect[2] = {u,v};
    double norm_sq = mn_euclid_norm(2,vect)*mn_euclid_norm(2,vect);
    while(norm_sq > 1 || fabs(norm_sq) == 0){
        u = gsl_rng_uniform(generator)*2.0 - 1.0;
        v = gsl_rng_uniform(generator)*2.0 - 1.0;
        vect[0] = u;
        vect[1] = v;
        norm_sq = mn_euclid_norm(2,vect)*mn_euclid_norm(2,vect);
    }
    double m = sqrt(-2.0*log(norm_sq)/norm_sq);
    result.x = u*m;
    result.y = v*m;

    return result;
}

//Euler-Maruyama
void mn_euler_maruyama_step(double t, double delta_t, double y[], SDEfunc func, int dimension, gsl_rng* gen, void *params){
    double *g = malloc(dimension*dimension*sizeof(double));
    double *W = malloc(dimension*sizeof(double));
    double *f = malloc(dimension*sizeof(double));
   
    func(t,y,f,g,params);

    //Zufallsvektor W erzeugen mit gaussian_random
    //stdabeichung sqrt(dt)
    for(int i = 0; i<dimension; i++){
        Tuple rnd_numb = mn_random_gaussian(gen);
        W[i] = (rnd_numb.x)*sqrt(delta_t);
    }
    //euler maruyama step
    for(int i = 0; i < dimension; i++){
        y[i] = y[i] + f[i]*delta_t + g[i]*W[i];
    }

    free(f);
    free(W);
    free(g);
}

/*--------------------------- Histogramme ------------------------------------*/

Hist* hist_malloc(int size, double low, double high){
    //Speicher allokieren
    Hist* H = malloc(sizeof(Hist));
    //Pfeile, weil nur pointer auf H, nicht H selber
    H->size = size;
    H->low = low;
    H->high = high;
    //Speicher für bin Array
    H->bin = malloc(size*sizeof(double));
    return H;
}

void hist_free(Hist* H){
    free(H->bin);
    free(H);
}

//Histogramm auf 0 setzen
void hist_reset(Hist* H){
    for(int i = 0; i < H->size; i++){
        H->bin[i] = 0.0;
    }
}

//Wert hinzufügen
void hist_add(Hist* H, double x){
    //Zahlen außerhalb von Grenzen werden ignoriert
    //Wenn x = H->high, dann liefer die erste formel bin
    //Nummer n+1, daher werden die zahlen in bin n gesteckt
    if(x < H->high && x >= H->low){
        int bin_number = floor(H->size*(x - H->low)/(H->high - H->low));
        H->bin[bin_number] += 1;
    }
    else if(x == H->high){
        H->bin[H->size - 1] += 1;
    }
}

//Ausgabe in CSV Datei
void hist_fprintf(Hist* H, FILE* file){
    //Breite von bin, ang. alle gleich
    double h = (H->high - H->low)/H->size;
    for(int i = 0; i < H->size; i++){
        //linke Grenze von bin, Dichte 
        fprintf(file, "%.16f, %.16f\n", H->low + i*h, H->bin[i]/h);
    }
}

//Anzahl von Samples
double hist_sum(Hist* H){
    double sum = 0.0;
    for(int i = 0; i < H->size; i++){
        sum += H->bin[i];
    }
    return sum;
}

//Ausgabe in CSV mit relativen Dichten
void hist_fprintf_rel(Hist* H, FILE* file){
    //Breite von bin, ang. alle gleich
    double h = (H->high - H->low)/(H->size);
    double N = hist_sum(H);
    for(int i = 0; i < H->size; i++){
        fprintf(file, "%.16f, %.16f\n", H->low + i*h, (H->bin[i]/h)/N);
    }
}

//Erwartungswert von Histogramm
double hist_mean(Hist* H){
    double N = hist_sum(H);
    //ang. alle bins gleich Groß, ist h die breite eines
    double h = (H->high - H->low)/H->size;
    //x_hat_0 ist Mittelpunkt des ersten bins
    double x_hat_0 = H->low + h/2.0;
    double sum = 0.0;
    for(int i = 0; i < H->size; i++){
        sum += (x_hat_0 + i*h)*(H->bin[i]);
    }
    return sum*(1.0/N);
}

//Standardabweichung von Histogramm
double hist_std(Hist* H){
    //Erwartungswert
    double E = hist_mean(H);
    //Summe aller Samples
    double N = hist_sum(H);
    //Breite von Bin
    double h = (H->high - H->low)/H->size;
    //Mittelpunkt von erten Bin
    double x_hat_0 = H->low + h/2.0;
    double sum = 0.0;
    for(int i = 0; i < H->size; i++){
        sum += pow((x_hat_0 + i*h - E),2.0)*(H->bin[i]);
    }
    sum = (1.0/(N - 1))*sum;
    return sqrt(sum);
}