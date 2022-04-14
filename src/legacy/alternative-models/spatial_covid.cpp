#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>

#include "solvers.h"
#include "real_data_reader.h"

#ifndef OUTPUT_PATH
#define OUTPUT_PATH "../output"
#endif

// Numeros de elementos de area
#ifndef NEA
#define NEA 10
#endif
// Numero de variaveis
#ifndef NA
#define NA 5
#endif

// Matriz F
static double F[NEA][NEA];

// Fluxos
// GamaS, GamaE, GamaY, GamaR, GamaN
static double GamaY_In[NEA][NA];
static double GamaY_Out[NEA][NA];

// Constante de proprocionalidade fluxo saida
static double csi = 0.0;

// Variaveis da EDO, valores ja divididos pro Ai
// si-ei-yi-ri-ni
//static double Y[NEA][NA];
//static double YOUT[NEA][NA];
static double Y_sum[NA];

// Area elemento de area
static double A = 10.0;
static double Ai[NEA];
// Area total
static double AT = A*NEA;

// dados 2019
double br_birth_rate = 13.82;
double br_death_rate = 6.503;
double br_population = 211049527.0;
double br_total_birth = br_birth_rate*br_population/1000.0;
double br_total_death = br_death_rate*br_population/1000.0;



// Taxa de Natalidade
//static double Lambda = 774835.0/(65640000.0*365.0);
static double Lambda = br_total_birth/(br_population*365);
static double lambda = Lambda/A;
// Taxa de contaminacao
static double beta = 3.4;
// Taxa de mortalidade equilibrio
//static double mu_eq = 774835.0/(65640000.0*365.0)/2.0;
static double mu_eq = br_total_death/(br_population * 365);
// Taxa de mortalidade covid
static double mu_cov = 0.01*mu_eq;
// Taxa de recuperacao
static double gama = 2.0/1924.0;
// Taxa de conversao de individuos infectados
static double a = 1.0/58.0;

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
    dydt[Var::S] = Lambda - mu_eq * y[Var::S] - beta * (y[Var::S]/y[Var::N]) * y[Var::Y];
    dydt[Var::E] = beta * (y[Var::S]/y[Var::N]) * y[Var::Y] - (mu_eq + a) * y[Var::E];
    dydt[Var::Y] = a * y[Var::E] - (mu_eq + mu_cov + gama) * y[Var::Y];
    dydt[Var::R] = gama * y[Var::Y] - mu_eq * y[Var::R];
    dydt[Var::N] = Lambda - mu_eq * y[Var::N] - mu_cov * y[Var::Y];
}

void derivs2D(double t, double y[][NA], double dydt[][NA])
{
    for (int i = 0; i < NEA; ++i)
    {
        dydt[i][Var::S] = (GamaY_In[i][Var::S] - GamaY_Out[i][Var::S]) / A + lambda - mu_eq * y[i][Var::S]
                        - beta * (y[i][Var::S]/y[i][Var::N]) * y[i][Var::Y];
        dydt[i][Var::E] = (GamaY_In[i][Var::E] - GamaY_Out[i][Var::E]) / A
                        + beta * (y[i][Var::S]/y[i][Var::N]) * y[i][Var::Y] - (mu_eq + a) * y[i][Var::E];
        dydt[i][Var::Y] = (GamaY_In[i][Var::Y] - GamaY_Out[i][Var::Y]) / A 
                        + a * y[i][Var::E] - (mu_eq + mu_cov + gama) * y[i][Var::Y];
        dydt[i][Var::R] = (GamaY_In[i][Var::R] - GamaY_Out[i][Var::R]) / A + gama * y[i][Var::Y] 
                        - mu_eq * y[i][Var::R];
        dydt[i][Var::N] = (GamaY_In[i][Var::N] - GamaY_Out[i][Var::N]) / A + lambda 
                        - mu_eq * y[i][Var::N] - mu_cov * y[i][Var::Y];
    }
}


// A soma será a da linha i da matrix f(Dij)
double sumF_n(int i) {
    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum += (i != j) ? F[i][j] * GamaY_Out[j][Var::N] : 0;
    }

    return sum;
}


double sumF_seir(double y[][NA], int i, Var::Enum type) {
    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum += (i != j) ? (y[j][type] / y[j][Var::N]) * F[i][j] * GamaY_Out[j][Var::N] : 0;
    }

    return sum;
}


void updateFlow(double y[][NA]) {
    for (int i = 0; i < NEA; ++i)
    {
        GamaY_Out[i][Var::N] = csi * y[i][Var::N];
        GamaY_Out[i][Var::S] = (y[i][Var::S] / y[i][Var::N]) * GamaY_Out[i][Var::N];
        GamaY_Out[i][Var::E] = (y[i][Var::E] / y[i][Var::N]) * GamaY_Out[i][Var::N];
        GamaY_Out[i][Var::Y] = (y[i][Var::Y] / y[i][Var::N]) * GamaY_Out[i][Var::N];
        GamaY_Out[i][Var::R] = (y[i][Var::R] / y[i][Var::N]) * GamaY_Out[i][Var::N];
    }

    for (int i = 0; i < NEA; ++i)
    {
        GamaY_In[i][Var::N]  = sumF_n(i);
        GamaY_In[i][Var::S]  = sumF_seir(y, i, Var::S);
        GamaY_In[i][Var::E]  = sumF_seir(y, i, Var::E);
        GamaY_In[i][Var::Y]  = sumF_seir(y, i, Var::Y);
        GamaY_In[i][Var::R]  = sumF_seir(y, i, Var::R);
    }
}

void driver_simple(double y0[], const char* filename)
{
    FILE *result = fopen(filename, "w");
	double h = 1.0, t=0.0, y[NA], dydt[NA], yout[NA];

    memcpy(y, y0, NA*sizeof(double));

    fprintf(result, "S,E,I,R,N\n");
    fprintf(result, "%f,%f,%f,%f,%f\n", y[Var::S], y[Var::E], y[Var::Y], y[Var::R], y[Var::N]);

    for (int i = 0; i < 800; i++)
    {
        derivs(t, y, dydt);
        rk4(y, dydt, t, h, yout, derivs);

        fprintf(result, "%f,%f,%f,%f,%f\n", yout[Var::S], yout[Var::E], yout[Var::Y], yout[Var::R], yout[Var::N]);
        t += h;
        memcpy(y, yout, NA*sizeof(double));
    }
    fclose(result);
}

void sumY(double y[][NA]) 
{
    Y_sum[Var::S] = 0;
    Y_sum[Var::E] = 0;
    Y_sum[Var::Y] = 0;
    Y_sum[Var::R] = 0;
    Y_sum[Var::N] = 0;

    for (int i = 0; i < NEA; ++i)
    {
        Y_sum[Var::S] += A * y[i][Var::S];
        Y_sum[Var::E] += A * y[i][Var::E];
        Y_sum[Var::Y] += A * y[i][Var::Y];
        Y_sum[Var::R] += A * y[i][Var::R];
        Y_sum[Var::N] += A * y[i][Var::N];
    }
}

void driver2D_simple(double y0[][NA], const char*filename)
{
    FILE *result = fopen(filename, "w");
	double h = 1.0, t=0.0, y[NEA][NA], dydt[NEA][NA], yout[NEA][NA];

    memcpy(y, y0, NEA*NA*sizeof(double));
    updateFlow(y);
    sumY(y);

    fprintf(result, "S,E,I,R,N\n");
    fprintf(result, "%f,%f,%f,%f,%f\n", Y_sum[Var::S], Y_sum[Var::E], Y_sum[Var::Y], Y_sum[Var::R], Y_sum[Var::N]);

    for (int i = 0; i < 700; i++)
    {
        derivs2D(t, y, dydt);
        rk42D(y, dydt, t, h, yout, derivs2D);

        memcpy(y, yout, NEA*NA*sizeof(double));
        updateFlow(y);

        sumY(y);
        fprintf(result, "%f,%f,%f,%f,%f\n", Y_sum[Var::S], Y_sum[Var::E], Y_sum[Var::Y], Y_sum[Var::R], Y_sum[Var::N]);
        t += h;
    }
    fclose(result);

}

int main()
{
    char file2open[4096];
	double h, t=0.0, y0[5];
    y0[Var::N] = br_population;
    y0[Var::E] = 0;
    y0[Var::Y] = 1;
    y0[Var::S] = y0[Var::N] - y0[Var::E] - y0[Var::Y];
    y0[Var::R] = 0;

    sprintf (file2open, "%s/result.txt", OUTPUT_PATH);
    driver_simple(y0, file2open);

    double Y0[NEA][NA];
    for (int i = 0; i < NEA; ++i)
    {
        Y0[i][Var::S] = y0[Var::S]/AT;
        Y0[i][Var::E] = 0;
        Y0[i][Var::Y] = 0;
        Y0[i][Var::R] = y0[Var::R]/AT;
        Y0[i][Var::N] = y0[Var::N]/AT;
        
        for (int j = 0; j < NEA; ++j)
        {
            //F[i][j] = (i != j) ? 1.0/(NEA * NEA - NEA) : 0.0;
            F[i][j] = 0;
        }
    }
    for (int j = 1; j < NEA; ++j)
    {
        F[0][j] =  1.0/(NEA - 1);
    }
    Y0[0][Var::Y] = y0[Var::Y];

    sprintf (file2open, "%s/result2D.txt", OUTPUT_PATH);
    driver2D_simple(Y0, file2open);


    return 0;
}
