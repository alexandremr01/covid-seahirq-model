//
// Created by Alexandre Maranhão on 14/04/22.
//

#ifndef COVID_SEAHIRQ_MODEL_PARAMETERS_H
#define COVID_SEAHIRQ_MODEL_PARAMETERS_H
#include "../seahirq_lib/macros.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <vector>

struct PhaseParameters {
    Eigen::Matrix<double, NEA, 1>  xI;
    Eigen::Matrix<double, NEA, 1>  xA;
    Eigen::Matrix<double, NEA, NEA> beta;
    PhaseParameters() {};
};

struct DynamicParameters {
public:
    DynamicParameters(
            Eigen::Matrix<double, NEA, 1>  xI0,
            Eigen::Matrix<double, NEA, 1>  xA0,
            Eigen::Matrix<double, NEA, NEA> beta0
    ){
        // std::cout << " Phase 0 " << phase0Parameters.xI;
        PhaseParameters *phase0Parameters = new PhaseParameters();
        phase0Parameters->xI = xI0;
        phase0Parameters->xA = xA0;
        phase0Parameters->beta = beta0;
        this->vectorPhaseParameters.push_back(phase0Parameters);
        this->days.push_back(0);
    };
    void addPhase(
            int startDay,
            Eigen::Matrix<double, NEA, 1>  xI,
            Eigen::Matrix<double, NEA, 1>  xA,
            Eigen::Matrix<double, NEA, NEA> beta
    ){
        // std::cout << " Phase 0 " << phase0Parameters.xI;
        PhaseParameters *phaseParameters = new PhaseParameters();
        phaseParameters->xI = xI;
        phaseParameters->xA = xA;
        phaseParameters->beta = beta;
        // std::cout << " Phase: " << phaseParameters->xI;
        this->vectorPhaseParameters.push_back(phaseParameters);
        this->days.push_back(startDay);
    };

    std::vector<PhaseParameters*> vectorPhaseParameters;
    std::vector<int> days;
};

struct StaticParameters {
    StaticParameters() {};

    Eigen::Matrix<double, NEA, 1> Lambda;
    // Taxa de mortalidade equilibrio
    Eigen::Matrix<double, NEA, 1> mu_eq;
    // Quanto a infectividade do assimptomatico é menor do que a do simptomatico(entre 0 e 1)
    Eigen::Matrix<double, NEA, 1> alpha;
    // Fator de correção para infecções pré-sintomáticas
    Eigen::Matrix<double, NEA, 1> ksi;
    // Probabilidade de um exposto se tornar sintomatico
    Eigen::Matrix<double, NEA, 1> rho;
    // Porcentagem de infectados que precisarao de hospitalizacao
    Eigen::Matrix<double, NEA, 1> phi;
    // Porcentagem de infectados que precisarao de leitos de UTI
    Eigen::Matrix<double, NEA, 1> eta;
    // Taxa de conversao de individuos infectados
    Eigen::Matrix<double, NEA, 1> a;
    // Taxa de mortalidade covid
    Eigen::Matrix<double, NEA, 1> mu_cov;
    // Coeficiente de taxa de recuperacao de assintomaticos (modelo SEAIR)
    Eigen::Matrix<double, NEA, 1> theta;
    // relação entre mortalidade de infectados e mortalidade de hospitalizados
    Eigen::Matrix<double, NEA, 1> gama_A;
    // Coeficiente de taxa de hospitalização de um infectado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_H;
    // Coeficiente de taxa de recuperação de hospitalizado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_HR;
    // Coeficiente de taxa de recuperação do infectado sintomático não quarentenado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RI;
    // Coeficiente de taxa de recuperação do infectado assintomático não quarentenado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RA;
    // Coeficiente de taxa de recuperação do quarentenado sintomático (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RQI;
    // Coeficiente de taxa de recuperação do quarentenado assintomático (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RQA;
    // Coeficiente de taxa de hospitalização do quarentenado sintomático (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_HQI;
    // Fração média de óbitos entre infectados
    Eigen::Matrix<double, NEA, 1> Tc;
    // Fração média de óbitos entre hospitalizados
    Eigen::Matrix<double, NEA, 1> Tlc;


    Eigen::Matrix<double, NEA, 1> gama;

    Eigen::Matrix<double, NEA, 1> rel_pop;
};

#endif //COVID_SEAHIRQ_MODEL_PARAMETERS_H
