#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "MyNumerics/my_numerics.h"

/*-------------------------------- Parameter ----------------------------------*/
const double S_0 = 70.0; //Aktienkurs in EUR
const double K = 100.0; //Ausübungspreis in EUR
const double T = 1.0; //Jahre
const double r = 0.12; //Risikofreier Zinssatz
const double sigm = 0.10; //Volatilität

/*---------------------------------- GSL RNG Set-Up ------------------------------*/
gsl_rng* get_gsl_gen(){
    gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(generator, time(NULL));
    return generator;
}

void free_gsl(gsl_rng* generator){
    if(generator != NULL){
        gsl_rng_free(generator);
    }
}

/*----------------------------- Auxilliary Functions ----------------------------*/
//Minimum von Array bestimmen
double array_min(double arr[], int length){
    //Minimum bestimmen
    double min = arr[0];
    for(int i = 1; i < length; i++){
        if(arr[i] <= min){
            min = arr[i];
        }
    }
    return min;
}

//Maximum von Array bestimmen
double array_max(double arr[], int length){
    double max = arr[0];
    for(int i = 1; i < length; i++){
        if(arr[i] >= max){
            max = arr[i];
        }
    }
    return max;
}

//geometrische Brownian motion Funktion
void gbm(gsl_rng* generator, double delta_t, double S[]){
    //Zwei Normalverteilte Zufallszahlen mit stdabw. 1, erzeugt mit Polarmethode
    //Nur eine benötigt, die andere wird verworfen
    Tuple rand_gaussian = mn_random_gaussian(generator);
    double Z = rand_gaussian.x;
    double arg = (r - (sigm*sigm)/2.0)*delta_t + sigm*sqrt(delta_t)*Z;
    S[0] = S[0]*exp(arg);
}

//d_1 und d_2 Parameter
double d_1(double t_init, double S_init){
    double factor = 1.0/(sigm*sqrt(T-t_init));
    double bracket = log(S_init/K) + (r + sigm*sigm/2.0)*(T-t_init);
    return factor*bracket;
}

double d_2(double t, double d1){
    return d1 - sigm*sqrt(T-t);
}

//CDF von Standard-Normalverteilung
double cdf_std_normal_distribution(double d){
    //benutzt die beziehung:
    // N(x) = 0.5*(1 + erf(x/sqrt(2))) 
    return 0.5*(1.0 + mn_erf_simpson(d/sqrt(2),1e-4));
}

//Call-Option Preis zum Zeitpunkt t
double call_price_analytic(double t_init, double S_init){
    double d1 = d_1(t_init, S_init);
    double N_d1 = cdf_std_normal_distribution(d1);
    double N_d2 = cdf_std_normal_distribution(d_2(t_init, d1));
    return fmax(S_init*N_d1 - K*exp(-r*(T-t_init))*N_d2,0.0);
}

/*-------------------------------- Monte-Carlo --------------------------------*/
double mc_step(gsl_rng* generator, double S[], int N){
    //Aktienverlauf berechnen und Preis zum Verfall zurückgeben
    double t = 0.0;
    double dt = 1.0/N;
    while(t<=T){
        gbm(generator,dt,S);
        t+=dt;
    }
    //Auszahlungsbetrag Barwert
    //die skalierung erfolg nach dem aufsummieren in mc_simulation()
    //double payoff = fmax(S[0] - K, 0.0);
    return S[0];
}

void mc_simulation(gsl_rng* generator, int N, double t_init, double S_init, double S[], FILE* data){
    //Simulationsmethode, hauptsächlich für Teil 4
    double callprice = 0.0;
    double var = 0.0;
    
    //Speicher für berechnung von Fehler allokieren
    double *Cp = malloc(N*sizeof(double));
    
    //Simulationsschleife
    for(int i = 0; i<N; i++){
        //Startwert der Aktie
        S[0] = S_init;
        //Zeitentwicklung der Aktie simulieren
        double S_T = mc_step(generator,S,N);
        double payoff = fmax(S_T - K, 0.0);
        callprice += payoff;
        Cp[i] = payoff*exp(-r*T);
    }

    //Mittlewert berechnen und skalieren
    callprice = (1.0/N)*callprice*exp(-r*T);

    //Fehler berechnen
    for(int j = 0; j<N; j++){
        var += (Cp[j] - callprice)*(Cp[j] - callprice);
    }
    var = sqrt(var/(N*(N-1)));

    //Ergebnisse in CSV Schreiben:
    //Aktienwert zum Zeitpunkt t_0, der Monte-Carlo Preis, der Fehler, die analytische Lsg
    fprintf(data, "%.16f, %.16f, %.16f, %.16f\n", S_init, callprice, var, call_price_analytic(t_init,S_init));

    free(Cp);
}

/*----------------------------------- Main ------------------------------------------*/

int main(){
    FILE* aktienkurs_dat = fopen("aktienkurs_data.csv", "w");
    FILE* analytic_price_time_data = fopen("analytic_call_prioce_time.csv", "w");
    FILE* ana_price_data = fopen("analytic_solution.csv", "w");
    FILE* pdf_data = fopen("aktien_pdf.csv", "w");
    FILE* N100_data = fopen("N100_data.csv", "w");
    FILE* N500_data = fopen("N500_data.csv", "w");
    FILE* N1000_data = fopen("N1000_data.csv", "w");
    FILE* N10000_data = fopen("N10000_data.csv", "w");

    double S[1] = {S_0};
    double dt = 0.01;
    double t = 0.0;
    double t_0 = 0.0;

    /*--------------------------------- TEIL 1 --------------------------------*/

    //GSL genetrator
    gsl_rng* gen = get_gsl_gen();

    //GBM Simulation
    while(t <= T){
        //den analytischen Preis berechnen
        double ana_call_price = call_price_analytic(t_0, S_0);
        //Daten in CSV eintragen
        fprintf(aktienkurs_dat, "%f, %.16f, %.16f\n", t, S[0], ana_call_price);
        //Den Wert von S_{t+dt} berechnen
        gbm(gen,dt,S);
        t += dt;
    }

    printf("Teil 1 fertig.\n");

    /*-------------------------------- TEIL 2 ---------------------------------*/
    int N = 10000;
    
    double callprice = 0.0;
    double var = 0.0;
    double *Cp = malloc(N*sizeof(double)); //für berechnung von Fehhler
    double *S_verfall = malloc(N*sizeof(double)); //Für Histogramm, raw S(T) Werte
    //Monte-Carlo Simulation mit N Schritten
    for(int i = 0; i<N; i++){
        //reset S
        S[0] = S_0;
        double S_T = mc_step(gen,S,N);
        S_verfall[i] = S_T;
        double payoff = fmax(S_T - K, 0.0);
        callprice += payoff*exp(-r*T);
        Cp[i] = payoff*exp(-r*T);
    }

    printf("Histogramm fuer Teil 3 generieren.\n");

    //Histogramm für PDF plot für Teil 3
    Hist* hist = hist_malloc(1000, array_min(S_verfall,N), array_max(S_verfall,N) + 1.0);
    hist_reset(hist);
    for(int i = 0; i<N; i++){
        hist_add(hist, S_verfall[i]);
    }
    hist_fprintf_rel(hist, pdf_data);
    printf("Erwartungswert von dem Histogramm ist %.5f EUR\n", hist_mean(hist));
    printf("Standardabweichung von dem Histogramm ist %.5f EUR\n", hist_std(hist));
    hist_free(hist);

    printf("Teil 3 fertig.\n");

    //callprice berechnen
    callprice = (1.0/N)*callprice;
    for(int j = 0; j<N; j++){
        var += (Cp[j] - callprice)*(Cp[j] - callprice);
    }
    var = sqrt(var/(N*(N-1)));

    free(Cp);
    free(S_verfall);

    printf("Monte-Carlo Callprice ist %.10f EUR mit Fehler %.10f EUR\n", callprice,var);
    printf("Analytischer Callprice ist %.10f EUR\n", call_price_analytic(t_0, S_0));
    printf("Teil 2 fertig.\n");

    /*----------------------------- TEIL 4 ----------------------------------*/
    double S4 = 70.0;
    double delta_S = 0.5;
    double S_end = 130.0;

    
    /*----------------- N = 100 ---------------*/
    N = 100;

    while(S4 <= S_end){
        mc_simulation(gen,N,t_0,S4,S,N100_data);
        S4 += delta_S;
    } 
    printf("N = %d fertig \n", N);

    /*----------------- N = 500 ---------------*/
    S4 = 70.0;
    N = 500;
    while(S4 <= S_end){
        mc_simulation(gen,N,t_0,S4,S,N500_data);
        fprintf(ana_price_data, "%.16f, %.16f\n", S4, call_price_analytic(t_0,S4));
        S4 += delta_S;
    } 
    printf("N = %d fertig \n", N);

    /*----------------- N = 1000 ---------------*/
    S4 = 70.0;
    N = 1000;
    while(S4 <= S_end){
        mc_simulation(gen,N,t_0,S4,S,N1000_data);
        S4 += delta_S;
    }
    printf("N = %d fertig \n", N);

    /*---------------- N = 10 000 ---------------*/
    S4 = 70.0;
    N = 10000;
    double avg_time = 0.0;
    while(S4<=S_end){
        clock_t beg,end;
        beg = clock();
        mc_simulation(gen,N,t_0,S4,S,N10000_data);
        end = clock();
        avg_time += (double) (end- beg)/CLOCKS_PER_SEC;
        printf("N 10 000: S = %f done\n", S4);
        S4 += delta_S;
    }
    avg_time = avg_time/120.0;

    printf("Durchschnittliche Zeit pro Startwert bei 10 000 Stuetzstellen ist %.2f Sekunden\n", avg_time);
    
    printf("Teil 4 fertig.\n");

    free_gsl(gen);
    fclose(ana_price_data);
    fclose(N100_data);
    fclose(N500_data);
    fclose(N1000_data);
    fclose(N10000_data);
    fclose(pdf_data);
    fclose(aktienkurs_dat);
    fclose(analytic_price_time_data);
    return 0;
}