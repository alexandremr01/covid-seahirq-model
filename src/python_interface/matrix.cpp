#include "matrix.h"
#include <iostream>

void read_array(double arr[NEA], Eigen::Matrix<double, NEA, 1> pyArr){
    for (int i = 0; i < NEA; i++) {
        arr[i] = pyArr(i, 0);
    }
}

void read_array_per_day(double array[MAX_DAYS][NEA], Eigen::Matrix<double, NEA, 1>  pyMatrix, int day) {
    for (int j = 0; j < NEA; j++)
        for (int k=day; k < MAX_DAYS; k++){
            array[k][j] =  pyMatrix(j, 0);
            // std::cout << "Preenchendo dia " << k << " faixa " << j << " com " << pyMatrix(j, 0);
        }
}

void read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], Eigen::Matrix<double, NEA, NEA>  pyMatrix, int day) {
    if (pyMatrix.rows() != NEA ||pyMatrix.rows() != NEA)
        throw std::exception();
    for (int i = 0; i < NEA; i++)
        for (int j = 0; j < NEA; j++)
            for (int k=day; k < MAX_DAYS; k++)
                matrix[k][i][j] =  pyMatrix(i, j);
}

void view_array(double arr[NEA], std::string name){
    std::cout<< name << ": ";
    for (int i = 0; i < NEA; i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}