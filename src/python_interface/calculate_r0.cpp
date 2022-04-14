#include "calculate_r0.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

double spectral_radius(Eigen::Matrix<double, 6*NEA, 6*NEA> matrix){
    double largest = 0;
    Eigen::EigenSolver<Eigen::Matrix<double, 6*NEA, 6*NEA>> solver(6*NEA);
    solver.compute(matrix, false);

    Eigen::VectorXcd eivals = solver.eigenvalues();
    int r = eivals.rows();
    for (int i=0; i<r; i++){
        std::complex<double> z = eivals(i);
        if (std::abs(z) > std::abs(largest))
            largest = std::abs(z);
    }
    return largest;
}

double calculateR0(StaticParameters *p, DynamicParameters *dp) {
    // # A unica entrada de novas infecções é em E, vindo de A, I ou E
    Eigen::Matrix<double, 6*NEA, 6*NEA>  F;
    F.setZero(6*NEA, 6*NEA);

    struct Mat {
        enum Enum
        {
            E, A, I, Qi, Qa, H
        };
    };
    for (int i=0; i < NEA; i++){
        for (int j=0; j < NEA; j++){
            double val = dp->vectorPhaseParameters.at(0)->beta(i, j) * p->rel_pop(i);
            F(Mat::E*NEA + i, j) = p->ksi(j) * val;
            F(Mat::E*NEA +i, NEA + j) = p->alpha(j) * val; // A
            F(Mat::E*NEA +i, 2*NEA + j) = val;
        }
    }
    // # Transições entre compartimentos de infectados
    Eigen::Matrix<double, 6*NEA, 6*NEA>  V;
    V.setZero(6*NEA, 6*NEA);

    double gama_QI[NEA];
    double gama_QA[NEA];

    for (int k = 0; k < NEA; k++)
    {
        gama_QI[k] = dp->vectorPhaseParameters.at(0)->xI[k] * (p->gama_H[k] + p->gama_RI[k]) / (1 - dp->vectorPhaseParameters.at(0)->xI[k]);
        gama_QA[k] = dp->vectorPhaseParameters.at(0)->xA[k] * p->gama_RA[k] / (1 - dp->vectorPhaseParameters.at(0)->xA[k]);
    }

    for (int i = 0; i < NEA; i++) {
        int j;
        // # Exposto -> Exposto
        j = Mat::E + i;
        V(j, j) =  p->mu_eq(i) +  p->a(i);
        // # Assintomático
        j = Mat::A*NEA + i;
        V(j, j) = p->mu_eq(i) +  p->gama_RA(i) +  gama_QA[i];
        V(j, Mat::E*NEA+i) = - p->a(i) * (1- p->rho(i));
        // # Infectado sintomático
        j = Mat::I*NEA + i;
        V(j, j) = p->mu_eq(i) + p->theta(i)*p->mu_cov(i) + p->gama_H(i) + p->gama_RI(i) + gama_QI[i];
        V(j, Mat::E*NEA+i) = -p->a(i)*p->rho(i);
        // # Quarentenado sintomático
        j = Mat::Qi*NEA + i;
        V(j, j) = p->mu_eq(i) + p->gama_HQI(i) + p->gama_RQI(i);
        V(j, Mat::I*NEA+i) = -gama_QI[i];
        // # Quarentenado assintomático
        j = Mat::Qa*NEA + i;
        V(j, j) = p->mu_eq(i) + p->gama_RQA(i);
        V(j, Mat::A*NEA+i) = -gama_QA[i];
        // # Hospitalizado
        j = Mat::H*NEA + i;
        V(j, j) = p->mu_eq(i) + p->gama_HR(i) + p->mu_cov(i);
        V(j, Mat::I*NEA+i) = -p->gama_H(i);
        V(j, Mat::Qi*NEA+i) = -p->gama_HQI(i);
    }

    // std::cout << "writing to file" << std::endl;
    // std::ofstream file1("F.c.csv");
    // file1 << F.format(CSVFormat);
    // file1.close();
    // std::ofstream file2("V.c.csv");
    // file2 << V.format(CSVFormat);
    // file2.close();

    Eigen::Matrix<double, 6*NEA, 6*NEA> next_generation_matrix = F * V.inverse();
    double R0 = spectral_radius(next_generation_matrix);

    return R0;
}
