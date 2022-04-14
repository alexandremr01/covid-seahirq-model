//
// Created by Alexandre Maranhão on 14/04/22.
//

#ifndef COVID_SEAHIRQ_MODEL_MATRIX_H
#define COVID_SEAHIRQ_MODEL_MATRIX_H
#include "../seahirq_lib/macros.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

void read_array(double arr[NEA], Eigen::Matrix<double, NEA, 1> pyArr);

void read_array_per_day(double array[MAX_DAYS][NEA], Eigen::Matrix<double, NEA, 1>  pyMatrix, int day);

void read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], Eigen::Matrix<double, NEA, NEA>  pyMatrix, int day);

void view_array(double arr[NEA], std::string name);

#endif //COVID_SEAHIRQ_MODEL_MATRIX_H
