#pragma once

#include "macros.h"
#include <Eigen/LU>


struct Var {
    enum Enum
    {
        S,
        E,
        A,
        H,
        I,
        R,
		Ri,
        Qi,
        Qa,
        N,
        C,
		D,
        L
    };
};

struct ScenarioOutput {
double Y_sum[NA][MAX_DAYS];
double YOUT[NEA][NA][MAX_DAYS];
};

struct ScenarioParameters {

Model::Enum model;

int day;
// Taxa de Natalidade

// ---- Parâmetros que fixos no tempo -----


double Lambda[NEA];
// Taxa de mortalidade equilibrio
double mu_eq[NEA];
// Quanto a infectividade do assimptomatico é menor do que a do simptomatico(entre 0 e 1)
double alpha[NEA];
// Fator de correção para infecções pré-sintomáticas
double ksi[NEA];
// Probabilidade de um exposto se tornar sintomatico
double rho[NEA];
// Porcentagem de infectados que precisarao de hospitalizacao
double phi[NEA];
// Porcentagem de infectados que precisarao de leitos de UTI
double eta[NEA];
// Taxa de conversao de individuos infectados
double a[NEA];
// Taxa de mortalidade covid
double mu_cov[NEA];
// Coeficiente de taxa de recuperacao de assintomaticos (modelo SEAIR)
double theta[NEA];
// relação entre mortalidade de infectados e mortalidade de hospitalizados
double gama_A[NEA];
// Coeficiente de taxa de hospitalização de um infectado (modelo SEAHIR-Qia)
double gama_H[NEA];
// Coeficiente de taxa de recuperação de hospitalizado (modelo SEAHIR-Qia)
double gama_HR[NEA];
// Coeficiente de taxa de recuperação do infectado sintomático não quarentenado (modelo SEAHIR-Qia)
double gama_RI[NEA];
// Coeficiente de taxa de recuperação do infectado assintomático não quarentenado (modelo SEAHIR-Qia)
double gama_RA[NEA];
// Coeficiente de taxa de recuperação do quarentenado sintomático (modelo SEAHIR-Qia)
double gama_RQI[NEA];
// Coeficiente de taxa de recuperação do quarentenado assintomático (modelo SEAHIR-Qia)
double gama_RQA[NEA];
// Coeficiente de taxa de hospitalização do quarentenado sintomático (modelo SEAHIR-Qia)
double gama_HQI[NEA];
// Fração média de óbitos entre infectados
double Tc[NEA];
// Fração média de óbitos entre hospitalizados
double Tlc[NEA];

// ---- Parâmetros que variam no tempo -----

// Taxa de contaminacao
double beta[MAX_DAYS][NEA][NEA];
// Taxa de recuperacao
double gama[MAX_DAYS][NEA];
// Coeficiente de taxa de movimento do infectado sintomático para a quarentena (modelo SEAHIR-Qia)
double gama_QI[MAX_DAYS][NEA];
// Coeficiente de taxa de movimento do infectado assintomático para a quarentena (modelo SEAHIR-Qia)
double gama_QA[MAX_DAYS][NEA];
// Probabilidade de testar um infectado sintomatico (modelo SEAHIR-Qia)
double xI[MAX_DAYS][NEA];
// Probabilidade de testar um infectado assintomatico (modelo SEAHIR-Qia)
double xA[MAX_DAYS][NEA];

};

void sumY(double y[][NA], ScenarioParameters *p, ScenarioOutput *output);
void driver2D_simple(DerivFunc derivs, double y0[][NA], ScenarioParameters *params, ScenarioOutput *output);
void driver2D_eigen(DerivFunc derivs, Eigen::Matrix<double, NA, NEA> y0, ScenarioParameters *p, ScenarioOutput *o);
//void run_benchmark(DerivFunc derivs, double beta_scalar, ScenarioParameters *params);
int write_output(const char* filename, ScenarioOutput *o);
