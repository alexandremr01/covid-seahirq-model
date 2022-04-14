#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include <cstring>

// Numeros de faixas etarias
#ifndef NEA
#define NEA 16
#endif

// Numero de variaveis
#ifndef NA
#define NA 9
#endif

#ifndef MAX_DAYS
#define MAX_DAYS 200
#endif

#ifndef MAX_DELTA
#define MAX_DELTA 20
#endif

#include "solvers.h"
#include "real_data_reader.h"

// Variaveis da EDO, valores ja divididos pro Ai
// si-ei-yi-ri-ni
const int BUFFER_SIZE = 4096;
static double Y_sum[NA][MAX_DAYS];
static double YOUT[NEA][NA][MAX_DAYS];

// dados 2019
double br_birth_rate = 13.82;
double br_death_rate = 6.503;
double br_population = 211049527.0;
double br_total_birth = br_birth_rate*br_population/1000.0;
double br_total_death = br_death_rate*br_population/1000.0;

// Taxa de Natalidade
static double Lambda[NEA] = { br_total_birth/(br_population*365) };
// Taxa de mortalidade equilibrio
static double mu_eq[NEA] = { br_total_death/(br_population * 365) };
// Tempo tipico de recuperacao
static int dI[NEA] = { 3 };
// Tempo tipico de obito
static int dC[NEA] = { 3 };
// Tempo tipico de incubacao
static int dL[NEA] = { 3 };
// Tempo tipico de inicio de infeccao
static int dH[NEA] = { 3 };
// Tempo tipico de alta
static int dA[NEA] = { 3 };
// Tempo tipico de inicio de infeccao até hospitalizacao
static int dU[NEA] = { 3 };
// Tempo tipico de alta da UTI
static int dAU[NEA] = { 3 };
// Tempo tipico de obito na UTI
static int dO[NEA] = { 3 };
// Quanto a infectividade do assimptomatico é menor do que a do simptomatico(entre 0 e 1)
static double alpha[NEA] = { 0.7 };
// Probabilidade de um exposto se tornar sintomatico
static double rho[NEA] = { 0.5 };
// Porcentagem de infectados que precisarao de hospitalizacao
static double psi[NEA];
// Porcentagem de infectados que precisarao de leitos de UTI
static double eta[NEA];

// Taxa de contaminacao
static double beta[MAX_DAYS][NEA][NEA] = { { {3.4} } };
// Taxa de mortalidade covid
static double mu_cov[MAX_DAYS][NEA] = { {0.1*mu_eq[0]} };
// Taxa de recuperacao
static double gama[MAX_DAYS][NEA] = { {2.0/1924.0} };

static int day = 0;
struct Var {
    enum Enum
    {
        S,
        E,
        A,
        I,
        R,
        N,
        C,
        H,
        L
    };
};

// A soma será a da linha i da matrix f(Dij)
double sum_beta(int k, double y[][NA]) {
    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum +=  beta[day][k][j]*(y[k][Var::S]/Y_sum[Var::N][day - 1])*y[j][Var::I];
    }

    return sum;
}

void derivs_sir(double t, double y[][NA], double dydt[][NA])
{
    for (int k = 0; k < NEA; ++k)
    {
        YOUT[k][Var::S][day] = y[k][Var::S];
        YOUT[k][Var::I][day] = y[k][Var::I];
        YOUT[k][Var::R][day] = y[k][Var::R];
        YOUT[k][Var::N][day] = y[k][Var::N];
        YOUT[k][Var::A][day] = y[k][Var::A];
        YOUT[k][Var::C][day] = y[k][Var::C];
        YOUT[k][Var::H][day] = y[k][Var::H];
        YOUT[k][Var::L][day] = y[k][Var::L];

        dydt[k][Var::S] = Lambda[k]*y[k][Var::N] - mu_eq[k] * y[k][Var::S] - sum_beta(k, y);
        dydt[k][Var::I] = sum_beta(k, y) 
            - mu_eq[k] * y[k][Var::I]
            - mu_cov[day][k]*(day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0)
            - gama[day][k]*(day - dI[k] >= 0 ? YOUT[k][Var::I][day - dI[k]] : 0);
        dydt[k][Var::R] = gama[day][k] * (day - dI[k] >= 0 ? YOUT[k][Var::I][day - dI[k]] : 0) 
            - mu_eq[k] * y[k][Var::R];
        dydt[k][Var::N] = Lambda[k]*y[k][Var::N] - mu_eq[k] * y[k][Var::N]
            - mu_cov[day][k] * (day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0);

        dydt[k][Var::C] = mu_cov[day][k] * (day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0);
        dydt[k][Var::H] = psi[k] * (
                  (day - dH[k] >= 0 ? YOUT[k][Var::I][day - dH[k]] : 0)
                - (day - dA[k] >= 0 ? YOUT[k][Var::H][day - dA[k]] : 0)
                + (day - dAU[k] >= 0 ? YOUT[k][Var::H][day - dAU[k]] : 0) * eta[k]
                );
    }
}

void derivs_seir(double t, double y[][NA], double dydt[][NA])
{
    for (int k = 0; k < NEA; ++k)
    {
        YOUT[k][Var::S][day] = y[k][Var::S];
        YOUT[k][Var::E][day] = y[k][Var::E];
        YOUT[k][Var::I][day] = y[k][Var::I];
        YOUT[k][Var::R][day] = y[k][Var::R];
        YOUT[k][Var::N][day] = y[k][Var::N];
        YOUT[k][Var::A][day] = y[k][Var::A];
        YOUT[k][Var::C][day] = y[k][Var::C];
        YOUT[k][Var::H][day] = y[k][Var::H];
        YOUT[k][Var::L][day] = y[k][Var::L];

        dydt[k][Var::S] = Lambda[k]*y[k][Var::N] - mu_eq[k] * y[k][Var::S] - sum_beta(k, y);

        dydt[k][Var::E] = sum_beta(k, y) 
            - mu_eq[k] * y[k][Var::E]
            - (day - dL[k] >= 0 ? YOUT[k][Var::E][day - dL[k]] : 0);

        dydt[k][Var::I] = (day - dL[k] >= 0 ? YOUT[k][Var::E][day - dL[k]] : 0) 
            - mu_eq[k] * y[k][Var::I]
            - mu_cov[day][k]*(day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0)
            - gama[day][k]*(day - dI[k] >= 0 ? YOUT[k][Var::I][day - dI[k]] : 0);

        dydt[k][Var::R] = gama[day][k] * (day - dI[k] >= 0 ? YOUT[k][Var::I][day - dI[k]] : 0) 
            - mu_eq[k] * y[k][Var::R];

        dydt[k][Var::N] = Lambda[k]*y[k][Var::N] - mu_eq[k] * y[k][Var::N]
            - mu_cov[day][k] * (day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0);

        dydt[k][Var::C] = mu_cov[day][k] * (day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0);
        dydt[k][Var::H] = psi[k] * (
                  (day - dH[k] >= 0 ? YOUT[k][Var::I][day - dH[k]] : 0)
                - (day - dA[k] >= 0 ? YOUT[k][Var::H][day - dA[k]] : 0)
                + (day - dAU[k] >= 0 ? YOUT[k][Var::H][day - dAU[k]] : 0) * eta[k]
                );
    }
}

// A soma será a da linha i da matrix f(Dij)
double sum_alpha_beta(int k, double y[][NA]) {
    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum +=  alpha[j]*beta[day][k][j]*(y[k][Var::S]/Y_sum[Var::N][day - 1])*y[j][Var::I];
    }

    return sum;
}

// A soma será a da linha i da matrix f(Dij)
double sum_beta_asint(int k, double y[][NA]) {
    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum +=  beta[day][k][j]*(y[k][Var::S]/Y_sum[Var::N][day - 1])*y[j][Var::A];
    }

    return sum;
}

void derivs_seair(double t, double y[][NA], double dydt[][NA])
{
    for (int k = 0; k < NEA; ++k)
    {
        YOUT[k][Var::S][day] = y[k][Var::S];
        YOUT[k][Var::E][day] = y[k][Var::E];
        YOUT[k][Var::A][day] = y[k][Var::A];
        YOUT[k][Var::I][day] = y[k][Var::I];
        YOUT[k][Var::R][day] = y[k][Var::R];
        YOUT[k][Var::N][day] = y[k][Var::N];
        YOUT[k][Var::C][day] = y[k][Var::C];
        YOUT[k][Var::H][day] = y[k][Var::H];
        YOUT[k][Var::L][day] = y[k][Var::L];

        dydt[k][Var::S] = Lambda[k]*y[k][Var::N] 
            - mu_eq[k] * y[k][Var::S] 
            - sum_beta_asint(k, y)
            - sum_alpha_beta(k, y);
        dydt[k][Var::E] = sum_beta_asint(k, y) 
            + sum_alpha_beta(k, y)
            - mu_eq[k] * y[k][Var::E]
            - (day - dL[k] >= 0 ? YOUT[k][Var::E][day - dL[k]] : 0);
        dydt[k][Var::I] = rho[k]*(day - dL[k] >= 0 ? YOUT[k][Var::E][day - dL[k]] : 0) 
            - mu_eq[k] * y[k][Var::I]
            - mu_cov[day][k]*(day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0)
            - gama[day][k]*(day - dI[k] >= 0 ? YOUT[k][Var::I][day - dI[k]] : 0);
        dydt[k][Var::A] = (1 - rho[k])*(day - dL[k] >= 0 ? YOUT[k][Var::E][day - dL[k]] : 0) 
            - mu_eq[k] * y[k][Var::A]
            - (day - dI[k] >= 0 ? YOUT[k][Var::A][day - dI[k]] : 0);
        dydt[k][Var::R] = gama[day][k] * (day - dI[k] >= 0 ? YOUT[k][Var::I][day - dI[k]] : 0) 
            + (day - dI[k] >= 0 ? YOUT[k][Var::A][day - dI[k]] : 0)
            - mu_eq[k] * y[k][Var::R];
        dydt[k][Var::N] = Lambda[k]*y[k][Var::N] 
            - mu_eq[k] * y[k][Var::N]
            - mu_cov[day][k] * (day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0);

        dydt[k][Var::C] = mu_cov[day][k] * (day - dC[k] >= 0 ? YOUT[k][Var::I][day - dC[k]] : 0);
        dydt[k][Var::H] = psi[k] * (
                  (day - dH[k] >= 0 ? YOUT[k][Var::I][day - dH[k]] : 0)
                - (day - dA[k] >= 0 ? YOUT[k][Var::H][day - dA[k]] : 0)
                + (day - dAU[k] >= 0 ? YOUT[k][Var::H][day - dAU[k]] : 0) * eta[k]
                );
    }
}

void sumY(double y[][NA])
{
    Y_sum[Var::S][day] = 0;
    Y_sum[Var::E][day] = 0;
    Y_sum[Var::I][day] = 0;
    Y_sum[Var::R][day] = 0;
    Y_sum[Var::N][day] = 0;
    Y_sum[Var::A][day] = 0;
    Y_sum[Var::C][day] = 0;
    Y_sum[Var::H][day] = 0;
    Y_sum[Var::L][day] = 0;

    for (int i = 0; i < NEA; ++i)
    {
        Y_sum[Var::S][day] +=  y[i][Var::S];
        Y_sum[Var::E][day] +=  y[i][Var::E];
        Y_sum[Var::I][day] +=  y[i][Var::I];
        Y_sum[Var::R][day] +=  y[i][Var::R];
        Y_sum[Var::N][day] +=  y[i][Var::N];
        Y_sum[Var::A][day] +=  y[i][Var::A];
        Y_sum[Var::C][day] +=  y[i][Var::C];
        Y_sum[Var::H][day] +=  y[i][Var::H];
        Y_sum[Var::L][day] +=  y[i][Var::L];
    }
}

void driver2D_simple(double y0[][NA], const char* filename)
{
    day = 0;
    FILE *result = fopen(filename, "w");
	double h = 1.0, t=0.0, y[NEA][NA], dydt[NEA][NA];
    //double yout[NEA][NA];
    double yout[NEA][NA];

    for (int k = 0; k < NEA; k++)
    {
        y[k][Var::S] =  y0[k][Var::S];
        y[k][Var::E] =  y0[k][Var::E];
        y[k][Var::I] =  y0[k][Var::I];
        y[k][Var::R] =  y0[k][Var::R];
        y[k][Var::N] =  y0[k][Var::N];
        y[k][Var::A] =  y0[k][Var::A];
        y[k][Var::C] =  y0[k][Var::C];
        y[k][Var::H] =  y0[k][Var::H];
        y[k][Var::I] =  y0[k][Var::I];
    }
    sumY(y);

    fprintf(result, "S,E,I,R,N,A,C,H,L\n");
    fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", 
              Y_sum[Var::S][day], Y_sum[Var::E][day], Y_sum[Var::I][day], Y_sum[Var::R][day], Y_sum[Var::N][day]
            , Y_sum[Var::A][day], Y_sum[Var::C][day], Y_sum[Var::H][day], Y_sum[Var::L][day]);
    for (int k = 0; k < NEA; k++)
    {
        YOUT[k][Var::S][day] = y[k][Var::S];
        YOUT[k][Var::E][day] = y[k][Var::E];
        YOUT[k][Var::I][day] = y[k][Var::I];
        YOUT[k][Var::R][day] = y[k][Var::R];
        YOUT[k][Var::N][day] = y[k][Var::N];
        YOUT[k][Var::A][day] = y[k][Var::A];
        YOUT[k][Var::C][day] = 0;
        YOUT[k][Var::H][day] = 0;
        YOUT[k][Var::L][day] = 0;
        fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", y[k][Var::S], y[k][Var::E], y[k][Var::I]
                , y[k][Var::R], y[k][Var::N], y[k][Var::A], y[k][Var::C], y[k][Var::H], y[k][Var::L]);
    }
    fprintf(result, "\n");

    for (int i = 1; i < MAX_DAYS; i++)
    {
        day = i;
        derivs_seir(t, y, dydt);
        rk42D(y, dydt, t, h, yout, derivs_seir);

        //memcpy(y, yout, NEA*NA*sizeof(double));
        for (int k = 0; k < NEA; k++)
        {
            y[k][Var::S] =  yout[k][Var::S];
            y[k][Var::E] =  yout[k][Var::E];
            y[k][Var::I] =  yout[k][Var::I];
            y[k][Var::R] =  yout[k][Var::R];
            y[k][Var::N] =  yout[k][Var::N];
            y[k][Var::A] =  yout[k][Var::A];
            y[k][Var::C] =  yout[k][Var::C];
            y[k][Var::H] =  yout[k][Var::H];
            y[k][Var::I] =  yout[k][Var::I];
        }

        sumY(y);
        fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", 
                  Y_sum[Var::S][day], Y_sum[Var::E][day], Y_sum[Var::I][day], Y_sum[Var::R][day], Y_sum[Var::N][day]
                , Y_sum[Var::A][day], Y_sum[Var::C][day], Y_sum[Var::H][day], Y_sum[Var::L][day]);
        for (int k = 0; k < NEA; k++)
        {
            YOUT[k][Var::S][day] = y[k][Var::S];
            YOUT[k][Var::E][day] = y[k][Var::E];
            YOUT[k][Var::I][day] = y[k][Var::I];
            YOUT[k][Var::R][day] = y[k][Var::R];
            YOUT[k][Var::N][day] = y[k][Var::N];
            YOUT[k][Var::A][day] = y[k][Var::A];
            YOUT[k][Var::C][day] = y[k][Var::C];
            YOUT[k][Var::H][day] = y[k][Var::H];
            YOUT[k][Var::L][day] = y[k][Var::L];
            fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", y[k][Var::S], y[k][Var::E], y[k][Var::I]
                    , y[k][Var::R], y[k][Var::N], y[k][Var::A], y[k][Var::C], y[k][Var::H], y[k][Var::L]);
        }
        fprintf(result, "\n");
        t += h;
    }
    fclose(result);
}

bool read_non_comment_line(char* buffer, FILE* input) {
  while(fgets(buffer, BUFFER_SIZE, input)) {
    if (*buffer != '#' && *buffer != '\n')
      return true;
  }
  return false;
}
bool read_simple_array(double array[NEA], char* buffer, FILE* input) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars=0;
  int chars_old;
  double val;
  char title[128];
  sscanf(buffer, "%s%n", title, &chars);
  printf("title: %s\n", title);
  for(int i=0; i < NEA; ++i) {
    sscanf(buffer + chars, "%lf%n", &val, &chars_old);
    chars += chars_old;
    array[i] = val;
    printf("  %lf", val);
  }
  printf("\n");
  return true;
}
bool read_simple_array(int array[NEA], char* buffer, FILE* input) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars=0;
  int chars_old;
  int val;
  char title[128];
  sscanf(buffer, "%s%n", title, &chars);
  printf("title: %s\n", title);
  for(int i=0; i < NEA; ++i) {
    sscanf(buffer + chars, "%d%n", &val, &chars_old);
    chars += chars_old;
    array[i] = val;
    printf("  %d", val);
  }
  printf("\n");
  return true;
}
bool read_array_per_day(double array[MAX_DAYS][NEA], char* buffer, FILE* input, int day) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars = 0;
  int chars_old;
  double val;
  char title[128];
  sscanf(buffer, "%s%n", title, &chars);
  printf("title: %s\n", title);
  for(int i=0; i < NEA; ++i) {
    sscanf(buffer + chars, "%lf%n", &val, &chars_old);
    chars += chars_old;
    for (int k = day; k < MAX_DAYS; ++k)
      {
        int idx = day == 0 ? 0 : k - 1;
        array[k][i] = val  >= 0 ? val  : array[idx][i];
      }
    printf("  %lf", val);
  }
  printf("\n");
  return true;
}
bool read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], char* buffer, FILE* input, int day) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars = 0;
  int chars_old;
  double val;
  char title[128];
  sscanf(buffer, "%s%n", title, &chars);
  printf("title: %s\n", title);
  for(int j = 0; j < NEA; ++j) {
    if (j>0) {
      read_non_comment_line(buffer, input);
      chars = 0;
    }
    for(int i = 0; i < NEA; ++i) {
      sscanf(buffer + chars, "%lf%n", &val, &chars_old);
      printf("  %lf", val);
      chars += chars_old;
      for (int k = day; k < MAX_DAYS; ++k) {
        int idx = day == 0 ? 0 : k - 1;
        matrix[k][j][i] = val  >= 0 ? val  : matrix[idx][j][i];
      }
    }
    printf("\n");
  }
  return true;
}

int main(int argc, char** argv)
{
    FILE *input;
    char output_filename[256];
    if (argc > 1)
    {
        input = fopen(argv[1], "r");
    }
    else
    {
        input = fopen("input/new-tempate.txt", "r");
    }
    if (argc > 2)
    {
      strncpy(output_filename, argv[2], 255);
    }
    else
    {
      strcpy(output_filename, "output/result2D.txt");
    }
    double h, t=0.0, y0[NEA][NA];
    double population[NEA], exposed[NEA], infected[NEA], removed[NEA];

    char buffer[BUFFER_SIZE];
    // Read Initial Values
    read_simple_array(population, buffer, input);
    read_simple_array(exposed, buffer, input);
    read_simple_array(infected, buffer, input);
    read_simple_array(removed, buffer, input);
    read_simple_array(Lambda, buffer, input);
    read_simple_array(mu_eq, buffer, input);
    read_simple_array(alpha, buffer, input);
    read_simple_array(rho, buffer, input);
    read_simple_array(psi, buffer, input);
    read_simple_array(eta, buffer, input);
    read_simple_array(dI, buffer, input);
    read_simple_array(dC, buffer, input);
    read_simple_array(dL, buffer, input);
    read_simple_array(dH, buffer, input);
    read_simple_array(dA, buffer, input);
    read_simple_array(dU, buffer, input);
    read_simple_array(dAU, buffer, input);
    read_simple_array(dO, buffer, input);

    for (int i = 0; i < NEA; ++i)
    {
        y0[i][Var::N] = population[i];
        y0[i][Var::E] = exposed[i];
        y0[i][Var::I] = infected[i];
        y0[i][Var::R] = removed[i];
        y0[i][Var::S] = y0[i][Var::N] - y0[i][Var::E] - y0[i][Var::I];
        y0[i][Var::A] = 0;
        y0[i][Var::C] = 0;
        y0[i][Var::H] = 0;
        y0[i][Var::L] = 0;
    }

    char title[128];
    while (read_non_comment_line(buffer, input))
    {
        sscanf(buffer, "%s %d", title, &day); //day
        printf("%s: %d\n", title, day);
        read_array_per_day(mu_cov, buffer, input, day);
        read_array_per_day(gama, buffer, input, day);
        read_matrix_per_day(beta, buffer, input, day);
    }
    fclose(input);

    driver2D_simple(y0, output_filename);
}
