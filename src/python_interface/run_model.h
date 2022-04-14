//
// Created by Alexandre Maranhão on 14/04/22.
//

#ifndef COVID_SEAHIRQ_MODEL_RUN_MODEL_H
#define COVID_SEAHIRQ_MODEL_RUN_MODEL_H
#include <vector>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "parameters.h"
#include "matrix.h"
#include <iostream>
#include "../seahirq_lib/macros.h"
#include "../seahirq_lib/scenario.h"
#include "../seahirq_lib/models.h"
#include "../seahirq_lib/solvers.h"

struct Output {
    std::vector<Eigen::Matrix<double, NA, 1>> Y_sum;  // for each day: # individuals per disease stage
    std::vector<Eigen::Matrix<double, NEA, NA>> YOUT; // for each day: # invidiuals per age strata and disease stage

    Output(double Y_sum[NA][MAX_DAYS], double YOUT[NEA][NA][MAX_DAYS]) {
        for (int i=0; i<MAX_DAYS; i++){
            Eigen::Matrix<double, NA, 1> Y_sum_daily;
            Eigen::Matrix<double, NEA, NA> YOUT_daily;
            for (int j=0; j<NA; j++){
                Y_sum_daily(j, 0) = Y_sum[j][i];
                for (int k=0; k<NEA; k++){
                    YOUT_daily(k, j) = YOUT[k][j][i];
                }
            }
            this->Y_sum.emplace_back(Y_sum_daily);
            this->YOUT.emplace_back(YOUT_daily);
        }
    };
    Output(){};
};


Output run_model(int model, Eigen::Matrix<double, NA, NEA>  y0, StaticParameters *p, DynamicParameters *dp, bool verbose);
#endif //COVID_SEAHIRQ_MODEL_RUN_MODEL_H
