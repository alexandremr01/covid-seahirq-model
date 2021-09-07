#pragma once
#include <stdio.h>
#include "macros.h"

struct DataReal
{
    int size;
    double data[MAX_DAYS_REAL];
};

DataReal reader(const char* csv_file, const char* country, char separator);

bool read_non_comment_line(char* buffer, FILE* input);

bool read_simple_array(double array[NEA], char* buffer, FILE* input);

bool read_simple_array(int array[NEA], char* buffer, FILE* input);

bool read_array_per_day(double array[MAX_DAYS][NEA], char* buffer, FILE* input, int day);

bool read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], char* buffer, FILE* input, int day);



