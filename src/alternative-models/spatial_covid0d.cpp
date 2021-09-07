#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>


// Numeros de elementos de area
#ifndef NEA
#define NEA 10
#endif
// Numero de variaveis
#ifndef NA
#define NA 5
#endif

#ifndef MAX_DAYS
#define MAX_DAYS 30
#endif


#include "solvers.h"
#include "real_data_reader.h"


// Area elemento de area
static double A = 10.0;
// Area total
static double AT = A*NEA;

// dados 2019
double br_birth_rate = 13.82;
double br_death_rate = 6.503;
double br_population = 211049527.0;
double br_total_birth = br_birth_rate*br_population/1000.0;
double br_total_death = br_death_rate*br_population/1000.0;



//// Taxa de Natalidade
////static double Lambda = 774835.0/(65640000.0*365.0);
//static double Lambda = br_total_birth/(br_population*365);
//static double lambda = Lambda/A;
//// Taxa de contaminacao
//static double beta = 3.4;
//// Taxa de mortalidade equilibrio
////static double mu_eq = 774835.0/(65640000.0*365.0)/2.0;
//static double mu_eq = br_total_death/(br_population * 365);
//// Taxa de mortalidade covid
//static double mu_cov = 0.1*mu_eq;
//// Taxa de recuperacao
//static double gama = 2.0/1924.0;
//// Taxa de conversao de individuos infectados
//static double a = 1.0/58.0;

// Taxa de Natalidade
//static double Lambda = 774835.0/(65640000.0*365.0);
static double Lambda[MAX_DAYS] = { br_total_birth/(br_population*365) };
static double lambda = { Lambda[0]/A };
// Taxa de contaminacao
static double beta[MAX_DAYS] = { 3.4 };
// Taxa de mortalidade equilibrio
//static double mu_eq = 774835.0/(65640000.0*365.0)/2.0;
static double mu_eq[MAX_DAYS] = { br_total_death/(br_population * 365) };
// Taxa de mortalidade covid
static double mu_cov[MAX_DAYS] = { 0.1*mu_eq[0] };
// Taxa de recuperacao
static double gama[MAX_DAYS] = { 2.0/1924.0 };
// Taxa de conversao de individuos infectados
static double a[MAX_DAYS] = { 1.0/58.0 };

static int day = 0;



struct Var {
    enum Enum
    {
        S,
        E,
        Y,
        R,
        N
    };
};

void derivs(double t, double y[], double dydt[])
{
    dydt[Var::S] = Lambda[day] - mu_eq[day] * y[Var::S] - beta[day] * (y[Var::S]/y[Var::N]) * y[Var::Y];
    dydt[Var::E] = beta[day] * (y[Var::S]/y[Var::N]) * y[Var::Y] - (mu_eq[day] + a[day]) * y[Var::E];
    dydt[Var::Y] = a[day] * y[Var::E] - (mu_eq[day] + mu_cov[day] + gama[day]) * y[Var::Y];
    dydt[Var::R] = gama[day] * y[Var::Y] - mu_eq[day] * y[Var::R];
    dydt[Var::N] = Lambda[day] - mu_eq[day] * y[Var::N] - mu_cov[day] * y[Var::Y];
}

void derivs_sir(double t, double y[], double dydt[])
{
    dydt[Var::S] = Lambda[day]*y[Var::N] - mu_eq[day] * y[Var::S] - beta[day] * (y[Var::S]/y[Var::N]) * y[Var::Y];
    dydt[Var::Y] = beta[day] * (y[Var::S]/y[Var::N]) * y[Var::Y] - (mu_eq[day] +mu_cov[day] + gama[day]) * y[Var::Y];
    dydt[Var::R] = gama[day] * y[Var::Y] - mu_eq[day] * y[Var::R];
    dydt[Var::N] = Lambda[day]*y[Var::N] - mu_eq[day] * y[Var::N] - mu_cov[day] * y[Var::Y];
}


void driver_simple(double y0[], const char* filename)
{
    FILE *result = fopen(filename, "w");
	double h = 1.0, t=0.0, y[NA], dydt[NA], yout[NA];

    memcpy(y, y0, NA*sizeof(double));

    //fprintf(result, "S,E,I,R,N\n");
    fprintf(result, "S,I,R,N\n");
    fprintf(result, "%f,%f,%f,%f\n", y[Var::S], y[Var::Y], y[Var::R], y[Var::N]);
    //fprintf(result, "%f,%f,%f,%f,%f\n", y[Var::S], y[Var::E], y[Var::Y], y[Var::R], y[Var::N]);

    for (int i = 0; i < MAX_DAYS; i++)
    {
        day = i;
        derivs_sir(t, y, dydt);
        rk4(y, dydt, t, h, yout, derivs_sir);

        //fprintf(result, "%f,%f,%f,%f,%f\n", yout[Var::S], yout[Var::E], yout[Var::Y], yout[Var::R], yout[Var::N]);
        fprintf(result, "%f,%f,%f,%f\n", yout[Var::S], yout[Var::Y], yout[Var::R], yout[Var::N]);
        t += h;
        memcpy(y, yout, NA*sizeof(double));
    }
    fclose(result);
}

int main(int argc, char** argv)
{
    FILE *input;
    char *output_filename;
    if (argc > 1)
    {
        input = fopen(argv[1], "r");
    }
    else
    {
        input = fopen("input/input.txt", "r");
    }
    if (argc > 2)
    {
      output_filename = argv[2];
    }
    else
    {
        output_filename = "output/result.txt";
    }
    double h, t=0.0, y0[5];
    double population, exposed, infected, removed;

    char buffer[300];
    // Read Initial Values
    while (fgets(buffer, sizeof(buffer), input)) {
        if (*buffer == '#' || *buffer == '\n') continue; /* Ignora comentários */
        else break;
    }
    if (sscanf(buffer, "%lf %lf %lf %lf\n", &population, &exposed, &infected, &removed) != 4) {
        printf("Invalid parameters input\n");
    }

    // Read Parameters
    while (fgets(buffer, sizeof(buffer), input)) {
        if (*buffer == '#' || *buffer == '\n') continue; /* Ignora comentários */
        else break;
    }
    if (sscanf(buffer, "%lf %lf %lf %lf %lf %lf\n", &Lambda[0], &beta[0], &mu_eq[0], &mu_cov[0], &gama[0], &a[0]) != 6) {
        printf("Invalid parameters input\n");
    }
    for (int i = 0; i < MAX_DAYS; ++i)
    {
        Lambda[i] = Lambda[0];
        beta[i]   = beta[0];
        mu_eq[i]  = mu_eq[0];
        mu_cov[i] = mu_cov[0];
        gama[i]   = gama[0];
        a[i]      = a[0];
    }
    int d = 1;
    while (fgets(buffer, sizeof(buffer), input) != NULL) {
        if (*buffer == '#' || *buffer == '\n') continue; /* Ignora comentários */
        int chars;
        sscanf(buffer, "%d%n", &day, &chars);
        sscanf(buffer+chars, "%lf %lf %lf %lf %lf %lf\n", &Lambda[day], &beta[day], &mu_eq[day], &mu_cov[day], &gama[day], &a[day]);
        for (int i = day; i < MAX_DAYS; ++i)
        {
            Lambda[i] = Lambda[day] >= 0 ? Lambda[day] : Lambda[i-1];
            beta[i]   = beta[day] >= 0 ? beta[day] : beta[i-1];
            mu_eq[i]  = mu_eq[day] >= 0 ? mu_eq[day] : mu_eq[i-1];
            mu_cov[i] = mu_cov[day] >= 0 ? mu_cov[day] : mu_cov[i-1];
            gama[i]   = gama[day] >= 0 ? gama[day] : gama[i-1];
            a[i]      = a[day] >= 0 ? a[day] : a[i-1];
        }
        d++;
    }
    fclose(input);

    y0[Var::N] = population;
    y0[Var::E] = exposed;
    y0[Var::Y] = infected;
    y0[Var::R] = removed;
    y0[Var::S] = y0[Var::N] - y0[Var::E] - y0[Var::Y];

    driver_simple(y0, output_filename);
    //printf("Lambda: %lf, beta: %lf, mu_eq: %lf, mu_cov: %lf, gama: %lf, a: %lf\n", Lambda, beta, mu_eq, mu_cov, gama, a);
    //FILE *china = fopen("output/china.txt", "w");
    FILE *brazil = fopen("output/brazil.txt", "w");
    //DataReal result = reader("data/time_series_covid19_confirmed_global.csv", "China");
    DataReal result2 = reader("data/time_series_covid19_confirmed_global.csv", "Brazil", ',');
    int i;
    //for (i = 0; i < result2.size; ++i)
    //{
        //if (result2.data[i] > 0.00001) break;
    //}

    //memmove(&result2.data[0], &result2.data[i], (result2.size - i) * sizeof(double));
    //result2.size -= i;
    //for (int i = 0; i < result.size; ++i)
    //{
        //fprintf(china, "%lf\n", result.data[i]);
    //}
    for (int i = 0; i < result2.size; ++i)
    {
        fprintf(brazil, "%lf\n", result2.data[i]);
    }
    //fclose(china);
    fclose(brazil);

    return 0;
}
