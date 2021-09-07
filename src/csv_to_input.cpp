#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include <stdbool.h>

#ifndef NEA
#define NEA 16
#endif

#define NUM_FAIXAS NEA

int main(int argc, char** argv) {
    FILE *csv;
    FILE *output;
    char filename[512];
	const char *input_path = "input/";
	char out_generated_path[512];

    if (argc == 2)
    {
        sprintf(filename, "%scenarios/%s/initial.csv", input_path, argv[1]);
        sprintf(out_generated_path, "%scenarios/%s", input_path, argv[1]);
		if ((csv = fopen( filename, "r")) == NULL)
		{
			perror("Abertura de arquivo de entrada initial.csv falhou. \n");
			return 1;
		}
    }
    else if (argc == 1)
    {
		sprintf(filename, "%scenarios/default/initial.csv", input_path);
		if ((csv = fopen(filename, "r")) == NULL)
		{
			perror("Abertura de arquivo de entrada initial.csv falhou. \n");
			return 1;
		}
    }
    else {
        printf("Uso: \n");
        printf("\tcsv_to_input [nome do cen[ario] \n");
        printf("\n");
        printf("Exemplos: \n");
        printf("1) Usando o cenário \"Default\":\n");
        printf("\tcsv_to_input \n");
        printf("2) Usando um cenário chamado de \"HardLockdown\":\n");
        printf("\tcsv_to_input HardLockdown\n");
        printf("\n");
        printf("Cenário: \n");
        printf("Um cenário é uma pasta dentro de input/cenarios com o nome do cenário e contendo os arquivos .csv referentes aos parâmetros da simulação que serão utilizados");
        printf("\t\n");
        return 0;
    }

    sprintf(filename, "%s/generated-input.txt", out_generated_path);
    if ((output = fopen(filename, "w")) == NULL)
	{
		perror("Abertura de arquivo de saida generated-input.txt falhou. \n");
		return 1;
	}

    fprintf(output, "# ---------------------------------------------------------- \n");
    fprintf(output, "#            INITIAL VALUES PER AGE GROUP                    \n");
    fprintf(output, "# ---------------------------------------------------------- \n");
    fprintf(output, "\n");

    char line[4096];
    while (fgets(line, 4096, csv) != NULL)
    {
        char* tmp = strdup(line);

        const char* tok;
		char* posn;
        for (tok = strtok_r(line, ",", &posn); tok && *tok; tok = strtok_r(NULL, ",\n", &posn))
        {
            fprintf(output, "%-15s\t", tok);
        }
        fprintf(output, "\n");
        free(tmp);
    }
    fclose(csv);

	// Open PARAMETERS file
    if (argc == 2)
    {
        sprintf(filename, "%scenarios/%s/parameters.csv", input_path, argv[1]);
        printf("Lendo: %s\n", filename);
        if((csv = fopen(filename, "r")) == NULL)
		{
			perror("Abertura de arquivo de parametros paramenters.csv falhou. \n");
			return 1;
		}
    }
    else if(argc == 1)
    {
        sprintf(filename, "%scenarios/default/parameters.csv", input_path);
		printf("Lendo: %s\n", filename);
		if((csv = fopen(filename, "r")) == NULL)
		{
			perror("Abertura de arquivo de parametros parameters.csv falhou. \n");
			return 1;
		}
    }
    else {
        printf("Uso: \n");
        printf("csv_to_input [nome do cenário] \n");
        printf("\t ou \n");
        printf("csv_to_input \n");
        printf("\t (Neste caso usará o cenário default) \n");
    }

	fprintf(output, "\n\n");
	fprintf(output, "# ---------------------------------------------------------- \n");
	fprintf(output, "#               PARAMETERS (fixed in time)                   \n");
	fprintf(output, "# ---------------------------------------------------------- \n");
	fprintf(output, "\n");

	fgets(line, 4096, csv); // Ignorar header
	while (fgets(line, 4096, csv) != NULL)
	{
		char* tmp = strdup(line);

		const char* tok;
		char* posn;
		for (tok = strtok_r(line, ",", &posn); tok && *tok; tok = strtok_r(NULL, ",\n", &posn))
		{
			fprintf(output, "%-15s\t", tok);
		}
		fprintf(output, "\n");
		free(tmp);
	}
	fprintf(output, "\n");
	fclose(csv);

	// Print gama, gama_QI, gama_QA and beta by day
	FILE *csv_beta_gama;

    if (argc == 2)
    {
        sprintf(filename, "%scenarios/%s/beta_gama.csv", input_path,  argv[1]);
        printf("Lendo: %s\n", filename);
		if ((csv_beta_gama = fopen(filename, "r")) == NULL)
		{
			perror("Abertura de arquivo de parametros beta_gama.csv falhou. \n");
			return 1;
		}
    }
    else if (argc == 1)
    {
        sprintf(filename, "%scenarios/default/beta_gama.csv", input_path);
		printf("Lendo: %s\n", filename);
		if ((csv_beta_gama = fopen(filename, "r")) == NULL)
		{
			perror("Abertura de arquivo de parametros beta_gama.csv falhou. \n");
			return 1;
		}
		sprintf(filename, "%scenarios/default/beta_gama.csv", input_path);
    }
    else {
        printf("Uso: \n");
        printf("csv_to_input [nome do cenário] \n");
        printf("\t ou \n");
        printf("csv_to_input \n");
        printf("\t (Neste caso usará o cenário default) \n");
    }

	// Print gama, gama_QI, gama_QA and beta by day
    fgets(line, 4096, csv_beta_gama); // Ignorar header
	int day = -1;

    while (fgets(line, 4096, csv_beta_gama) != NULL)
    {
		char* beta_tmp = strdup(line);
		bool first_tok = true;
        const char* tok;
		char *posn;
		for (tok = strtok_r(line, ",", &posn); tok && *tok; tok = strtok_r(NULL, ",\n", &posn))
		{
			if (first_tok) {
				if (atoi(tok) > day) {
					day = atoi(tok);
					fprintf(output, "# ---------------------------------------------------------- \n");
					fprintf(output, "DAY %s \n", tok);

					// Print GAMA
					fprintf(output, "GAMA \t\t");
					tok = strtok_r(NULL, ",\n", &posn); // Step
					for (int ii = 0; ii < NEA; ii++) {
						fprintf(output, "%-18s ", tok);
						tok = strtok_r(NULL, ",\n", &posn); // Step
					}
					fprintf(output, "\n");

					// Print xI
					fprintf(output, "xI \t\t\t");
					for (int ii = 0; ii < NEA; ii++) {
						fprintf(output, "%-18s ", tok);
						tok = strtok_r(NULL, ",\n", &posn); // Step
					}
					fprintf(output, "\n");

					// Print xA
					fprintf(output, "xA \t\t\t");
					for (int ii = 0; ii < NEA; ii++) {
						fprintf(output, "%-18s ", tok);
						if (ii != (NEA - 1)) tok = strtok_r(NULL, ",\n", &posn); // Step
					}
					fprintf(output, "\n");

					// Print BETA
					fprintf(output, "BETA \t\t");
				}
				else if (tok[0] != 10) {
					// Padding
					fprintf(output, "     	    ");
				}
				first_tok = false;
			}
			else {
				fprintf(output, "%-18s ", tok);
			}	
		}
		first_tok = true;
		fprintf(output, "\n");
		free(beta_tmp);
    }
    fclose(csv_beta_gama);


    return 0;
}
