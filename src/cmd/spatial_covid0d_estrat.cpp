#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>

#include "../seahirq_lib/scenario.h"
#include "../seahirq_lib/models.h"

#include "../seahirq_lib/macros.h"
#include "../seahirq_lib/solvers.h"
#include "real_data_reader.h"

// Variaveis da EDO, valores ja divididos pro Ai
// si-ei-yi-ri-ni
const int BUFFER_SIZE = 4096;
int main(int argc, char** argv)
{
    FILE *input;
    char output_filename[416];
    Model::Enum model = Model::SEAIR;
    if (argc > 1)
    {
		if ((input = fopen(argv[1], "r")) == NULL) {
			perror("Abertura de arquivo falhou. \n");
			return 1;
		}
    }
    else
    {
        printf("Uso: \n");
        printf("\tspatial_covid0d_estrat input_file [output_file] \n");
        printf("\nExemplos: \n");
        printf("\tspatial_covid0d_estrat input\\dados.txt\n");
        printf("\tspatial_covid0d_estrat input\\dados.txt output\\resultado.txt\n");
        return 0;
    }
    if (argc > 2)
    {
      strncpy(output_filename, argv[2], 255);
    }
    else
    {
      sprintf(output_filename, "output\\result.csv");
    }
    if (argc > 3)
    {
        int m = atoi(argv[3]);
        if (0 <= m && m <= Model::SEAHIRQ)
            model = (Model::Enum)atoi(argv[3]); 
        else {
            printf("Modelo invalido! \n");
            return 0;
        }
    }
    double y0[NEA][NA];
	double population[NEA], exposed[NEA], asymptomatic[NEA], infected[NEA], removed[NEA], recovered_symptomatic[NEA];
    ScenarioParameters params;
    ScenarioOutput output;
    params.day = 0;

    char buffer[BUFFER_SIZE];
    // Read Initial Values
    read_simple_array(population, buffer, input);
    read_simple_array(exposed, buffer, input);
    read_simple_array(infected, buffer, input);
    read_simple_array(removed, buffer, input);
	read_simple_array(recovered_symptomatic, buffer, input);
	read_simple_array(asymptomatic, buffer, input);
    read_simple_array(params.Lambda, buffer, input);
    read_simple_array(params.mu_eq, buffer, input);
    read_simple_array(params.mu_cov, buffer, input);
	read_simple_array(params.theta, buffer, input);
    read_simple_array(params.alpha, buffer, input);
    read_simple_array(params.ksi, buffer, input);
    read_simple_array(params.rho, buffer, input);
    read_simple_array(params.phi, buffer, input);
    read_simple_array(params.eta, buffer, input);
    read_simple_array(params.a, buffer, input);
    read_simple_array(params.gama_H, buffer, input);
    read_simple_array(params.gama_HR, buffer, input);
    read_simple_array(params.gama_RI, buffer, input);
    read_simple_array(params.gama_RA, buffer, input);
    read_simple_array(params.gama_RQI, buffer, input);
    read_simple_array(params.gama_RQA, buffer, input);
    read_simple_array(params.gama_HQI, buffer, input);
	read_simple_array(params.Tc, buffer, input);
	read_simple_array(params.Tlc, buffer, input);

    for (int i = 0; i < NEA; ++i)
    {
        y0[i][Var::N] = population[i];
        y0[i][Var::E] = exposed[i];
        y0[i][Var::I] = infected[i];
        y0[i][Var::R] = removed[i];
		y0[i][Var::Ri] = recovered_symptomatic[i];
        y0[i][Var::S] = y0[i][Var::N] - y0[i][Var::E] - y0[i][Var::I];
        y0[i][Var::A] = 0;
        y0[i][Var::C] = 0;
		y0[i][Var::D] = 0;
		y0[i][Var::A] = asymptomatic[i];
        y0[i][Var::H] = 0;
        y0[i][Var::L] = 0;
        y0[i][Var::Qi] = 0;
        y0[i][Var::Qa] = 0;
    }

    char title[128];
    while (read_non_comment_line(buffer, input))
    {
		if (buffer[0] == ' ') continue;
		sscanf(buffer, "%s %d", title,  &params.day); //day
        read_array_per_day(params.gama, buffer, input, params.day);
		read_array_per_day(params.xI, buffer, input, params.day);
		read_array_per_day(params.xA, buffer, input, params.day);
		read_matrix_per_day(params.beta, buffer, input, params.day);
    }
    fclose(input);

	for (int k = 0; k < NEA; k++)
	{
		for (int day = 0; day < MAX_DAYS; day++)
		{
			params.gama_QI[day][k] = params.xI[day][k] * (params.gama_H[k] + params.gama_RI[k]) / (1 - params.xI[day][k]);
			params.gama_QA[day][k] = params.xA[day][k] * params.gama_RA[k] / (1 - params.xA[day][k]);
		}
	}
	
    params.model = model;
    DerivFunc derivs = get_model(params.model);
    driver2D_simple(derivs, y0, &params, &output);
	if (write_output(output_filename, &output) == 1) {
		return 1;
	}
	else {
		printf("Exececucao bem sucedida. Verificar consistencia dos resultados.\n");
		return 0;
	}
    //run_benchmark(derivs, 5.0, &params);
}