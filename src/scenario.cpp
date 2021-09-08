#include "scenario.h"
#include <Eigen/Dense>
#include <iostream>
#include "models.h"
#include "solvers.h"


using namespace Eigen;

void sumY(double y[][NA], ScenarioParameters *p, ScenarioOutput *o)
{
    int d = MAX_DAYS;
    double (*Y_sum)[MAX_DAYS] = o->Y_sum;

    Y_sum[Var::S][p->day] = 0;
    Y_sum[Var::E][p->day] = 0;
    Y_sum[Var::I][p->day] = 0;
    Y_sum[Var::R][p->day] = 0;
	Y_sum[Var::Ri][p->day] = 0;
    Y_sum[Var::N][p->day] = 0;
    Y_sum[Var::A][p->day] = 0;
    Y_sum[Var::C][p->day] = 0;
    Y_sum[Var::H][p->day] = 0;
	Y_sum[Var::D][p->day] = 0;
    Y_sum[Var::L][p->day] = 0;
    Y_sum[Var::Qi][p->day] = 0;
    Y_sum[Var::Qa][p->day] = 0;

    for (int i = 0; i < NEA; ++i)
    {
        Y_sum[Var::S][p->day] +=  y[i][Var::S];
        Y_sum[Var::E][p->day] +=  y[i][Var::E];
        //Y_sum[Var::I][p->day] +=  y[i][Var::I] + y[i][Var::A] + y[i][Var::R] + y[i][Var::C];
        Y_sum[Var::I][p->day] +=  y[i][Var::I];
        Y_sum[Var::R][p->day] +=  y[i][Var::R];
		Y_sum[Var::Ri][p->day] += y[i][Var::Ri];
        Y_sum[Var::N][p->day] +=  y[i][Var::N];
        Y_sum[Var::A][p->day] +=  y[i][Var::A];
        Y_sum[Var::C][p->day] +=  y[i][Var::C];
		Y_sum[Var::D][p->day] += y[i][Var::D];
        Y_sum[Var::H][p->day] +=  y[i][Var::H];
        Y_sum[Var::L][p->day] +=  y[i][Var::L];
        Y_sum[Var::Qi][p->day] +=  y[i][Var::Qi];
        Y_sum[Var::Qa][p->day] +=  y[i][Var::Qa];
    }
}


void driver2D_simple(DerivFunc derivs, double y0[][NA], ScenarioParameters *p, ScenarioOutput *o)
{
    const int STEPS_PER_DAY = 10;
    p->day = 0;
    double (*YOUT)[NA][MAX_DAYS] = o->YOUT;
    double (*Y_sum)[MAX_DAYS] = o->Y_sum;

    double h = 1, t=0.0, y[NEA][NA] = { { 0 } }, dydt[NEA][NA];
    double yout[NEA][NA];

    for (int k = 0; k < NEA; k++) {
        y[k][Var::S] =  y0[k][Var::S];
        y[k][Var::E] =  y0[k][Var::E];
        y[k][Var::A] =  y0[k][Var::A];
        y[k][Var::H] =  y0[k][Var::H];
		y[k][Var::D] = y0[k][Var::D];
		y[k][Var::I] =  y0[k][Var::I];
        y[k][Var::R] =  y0[k][Var::R];
		y[k][Var::Ri] = y0[k][Var::Ri];
        y[k][Var::Qi] =  y0[k][Var::Qi];
        y[k][Var::Qa] =  y0[k][Var::Qa];
        y[k][Var::N] =  y0[k][Var::N];
        y[k][Var::C] =  y0[k][Var::C];
        y[k][Var::L] =  y0[k][Var::L];
    }
    sumY(y, p, o);
	
    for (int k = 0; k < NEA; k++)
    {
        YOUT[k][Var::S][p->day] = y[k][Var::S];
        YOUT[k][Var::E][p->day] = y[k][Var::E];
        YOUT[k][Var::I][p->day] = y[k][Var::I];
        YOUT[k][Var::R][p->day] = y[k][Var::R];
		YOUT[k][Var::Ri][p->day] = y[k][Var::Ri];
        YOUT[k][Var::N][p->day] = y[k][Var::N];
        YOUT[k][Var::A][p->day] = y[k][Var::A];
        YOUT[k][Var::H][p->day] = y[k][Var::H];
		YOUT[k][Var::D][p->day] = y[k][Var::D];
        YOUT[k][Var::Qi][p->day] = y[k][Var::Qi];
        YOUT[k][Var::Qa][p->day] = y[k][Var::Qa];
        YOUT[k][Var::C][p->day] = y[k][Var::C];
        YOUT[k][Var::L][p->day] = y[k][Var::L];
    }

    for (int i = 1; i < MAX_DAYS; i++)
    {
        float t2=t;
        p->day = i;
        for (int i=0; i < STEPS_PER_DAY; i++){
            derivs(t2, y, dydt, p);
            rk42D(y, dydt, t2, h/STEPS_PER_DAY, yout, derivs, p);
            //memcpy(y, yout, NEA*NA*sizeof(double));
            for (int k = 0; k < NEA; k++)
            {
                y[k][Var::S] =  yout[k][Var::S];
                y[k][Var::E] =  yout[k][Var::E];
                y[k][Var::I] =  yout[k][Var::I];
                y[k][Var::R] =  yout[k][Var::R];
                y[k][Var::Ri] = yout[k][Var::Ri];
                y[k][Var::N] =  yout[k][Var::N];
                y[k][Var::A] =  yout[k][Var::A];
                y[k][Var::C] =  yout[k][Var::C];
                y[k][Var::D] = yout[k][Var::D];
                y[k][Var::H] =  yout[k][Var::H];
                y[k][Var::L] =  yout[k][Var::L];
                y[k][Var::Qi] = yout[k][Var::Qi];
                y[k][Var::Qa] = yout[k][Var::Qa];
            }
            update_params(y, p);
            t2 += h/STEPS_PER_DAY;
        }
        sumY(y, p, o);
        
        for (int k = 0; k < NEA; k++)
        {
            YOUT[k][Var::S][p->day] = y[k][Var::S];
            YOUT[k][Var::E][p->day] = y[k][Var::E];
            YOUT[k][Var::I][p->day] = y[k][Var::I];
            YOUT[k][Var::R][p->day] = y[k][Var::R];
			YOUT[k][Var::Ri][p->day] = y[k][Var::Ri];
            YOUT[k][Var::N][p->day] = y[k][Var::N];
            YOUT[k][Var::A][p->day] = y[k][Var::A];
            YOUT[k][Var::C][p->day] = y[k][Var::C];
            YOUT[k][Var::H][p->day] = y[k][Var::H];
			YOUT[k][Var::D][p->day] = y[k][Var::D];
            YOUT[k][Var::L][p->day] = y[k][Var::L];
            YOUT[k][Var::Qi][p->day] = y[k][Var::Qi];
            YOUT[k][Var::Qa][p->day] = y[k][Var::Qa];
        }
        t += h;
    }
}

void driver2D_eigen(DerivFunc derivs, Eigen::Matrix<double, NA, NEA> y0, ScenarioParameters *p, ScenarioOutput *o)
{
    const int STEPS_PER_DAY = 10;
    p->day = 0;
    double (*YOUT)[NA][MAX_DAYS] = o->YOUT;
    double (*Y_sum)[MAX_DAYS] = o->Y_sum;

    double h = 1, t=0.0, y[NEA][NA] = { { 0 } }, dydt[NEA][NA];
    double yout[NEA][NA];

    for (int k = 0; k < NEA; k++) {
        y[k][Var::S] =  y0(Var::S, k);
        y[k][Var::E] =  y0(Var::E, k);
        y[k][Var::A] =  y0(Var::A, k);
        y[k][Var::H] =  y0(Var::H, k);
		y[k][Var::D] = y0(Var::D, k);
		y[k][Var::I] =  y0(Var::I, k);
        y[k][Var::R] =  y0(Var::R, k);
		y[k][Var::Ri] = y0(Var::Ri, k);
        y[k][Var::Qi] =  y0(Var::Qi, k);
        y[k][Var::Qa] =  y0(Var::Qa, k);
        y[k][Var::N] =  y0(Var::N, k);
        y[k][Var::C] =  y0(Var::C, k);
        y[k][Var::L] =  y0(Var::L, k);
    }
    sumY(y, p, o);
	
    for (int k = 0; k < NEA; k++)
    {
        YOUT[k][Var::S][p->day] = y[k][Var::S];
        YOUT[k][Var::E][p->day] = y[k][Var::E];
        YOUT[k][Var::I][p->day] = y[k][Var::I];
        YOUT[k][Var::R][p->day] = y[k][Var::R];
		YOUT[k][Var::Ri][p->day] = y[k][Var::Ri];
        YOUT[k][Var::N][p->day] = y[k][Var::N];
        YOUT[k][Var::A][p->day] = y[k][Var::A];
        YOUT[k][Var::H][p->day] = y[k][Var::H];
		YOUT[k][Var::D][p->day] = y[k][Var::D];
        YOUT[k][Var::Qi][p->day] = y[k][Var::Qi];
        YOUT[k][Var::Qa][p->day] = y[k][Var::Qa];
        YOUT[k][Var::C][p->day] = y[k][Var::C];
        YOUT[k][Var::L][p->day] = y[k][Var::L];
    }

    for (int i = 1; i < MAX_DAYS; i++)
    {
        float t2=t;
        p->day = i;
        for (int i=0; i < STEPS_PER_DAY; i++){
            derivs(t2, y, dydt, p);
            rk42D(y, dydt, t2, h/STEPS_PER_DAY, yout, derivs, p);
            //memcpy(y, yout, NEA*NA*sizeof(double));
            for (int k = 0; k < NEA; k++)
            {
                y[k][Var::S] =  yout[k][Var::S];
                y[k][Var::E] =  yout[k][Var::E];
                y[k][Var::I] =  yout[k][Var::I];
                y[k][Var::R] =  yout[k][Var::R];
                y[k][Var::Ri] = yout[k][Var::Ri];
                y[k][Var::N] =  yout[k][Var::N];
                y[k][Var::A] =  yout[k][Var::A];
                y[k][Var::C] =  yout[k][Var::C];
                y[k][Var::D] = yout[k][Var::D];
                y[k][Var::H] =  yout[k][Var::H];
                y[k][Var::L] =  yout[k][Var::L];
                y[k][Var::Qi] = yout[k][Var::Qi];
                y[k][Var::Qa] = yout[k][Var::Qa];
            }
            update_params(y, p);
            t2 += h/STEPS_PER_DAY;
        }
        sumY(y, p, o);
        
        for (int k = 0; k < NEA; k++)
        {
            YOUT[k][Var::S][p->day] = y[k][Var::S];
            YOUT[k][Var::E][p->day] = y[k][Var::E];
            YOUT[k][Var::I][p->day] = y[k][Var::I];
            YOUT[k][Var::R][p->day] = y[k][Var::R];
			YOUT[k][Var::Ri][p->day] = y[k][Var::Ri];
            YOUT[k][Var::N][p->day] = y[k][Var::N];
            YOUT[k][Var::A][p->day] = y[k][Var::A];
            YOUT[k][Var::C][p->day] = y[k][Var::C];
            YOUT[k][Var::H][p->day] = y[k][Var::H];
			YOUT[k][Var::D][p->day] = y[k][Var::D];
            YOUT[k][Var::L][p->day] = y[k][Var::L];
            YOUT[k][Var::Qi][p->day] = y[k][Var::Qi];
            YOUT[k][Var::Qa][p->day] = y[k][Var::Qa];
        }
        t += h;
    }
}

int cmp(const void *a, const void *b)
{
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
 
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

/*void run_trial(DerivFunc derivs, ScenarioParameters *params, double dI, double I0, double min_r0, 
                double max_r0, Matrix<double, NEA, NEA, RowMajor> C_esc, const char* output_filename)
{
    double t=0.0, y0[NEA][NA];
    ScenarioOutput output;
    double gama_escalar = 1 - exp(-1.0/dI);
    double R0 = min_r0;
    double step = (max_r0 - min_r0)/200.0;

    double bench_out[MAX_DAYS][200];
    for (int i = 0; i < 200; i++)
    {
        params->day = 0;
        for (int ii = 0; ii < NEA; ++ii)
        {
            y0[ii][Var::N] = p_vec[ii];
            y0[ii][Var::E] = 0;
            y0[ii][Var::I] = 0;
            y0[ii][Var::R] = 0;
			y0[ii][Var::Ri] = 0;
            y0[ii][Var::S] = y0[ii][Var::N] - y0[ii][Var::E] - y0[ii][Var::I];
            y0[ii][Var::A] = 0;
            y0[ii][Var::C] = 0;
            y0[ii][Var::H] = 0;
            y0[ii][Var::L] = 0;
            y0[ii][Var::Qi] = 0;
            y0[ii][Var::Qa] = 0;
        }
        y0[NEA/2][Var::I] = I0;
        y0[NEA/2][Var::S] = y0[NEA/2][Var::N] - y0[NEA/2][Var::E] - y0[NEA/2][Var::I];
        double beta_scalar = gama_escalar * R0;
        Matrix<double, NEA, NEA, RowMajor> b = beta_scalar * C_esc;
        for (int d = 0; d < MAX_DAYS; d++) 
        {
            for (int ii = 0; ii < NEA; ii++)
            {
                params->gama[d][ii] = gama_escalar;
                params->gama_A[ii] = gama_escalar;
                for (int jj = 0; jj < NEA; jj++)
                {
                    params->beta[d][ii][jj] = b(ii, jj);
                }
            }
        }
        driver2D_simple(derivs, y0, params, &output);
        for (int d = 0; d < MAX_DAYS; d++)
        {
            bench_out[d][i] = output.Y_sum[Var::I][d];
            bench_out[d][i] += output.Y_sum[Var::A][d];
            bench_out[d][i] += output.Y_sum[Var::R][d];
            bench_out[d][i] += output.Y_sum[Var::C][d];
        }

        R0 += i*step;
    }
    for (int d = 0; d < MAX_DAYS; d++)
    {
        qsort(bench_out[d], 200, sizeof(double), cmp);
        output.Y_sum[Var::I][d] = (bench_out[d][99] + bench_out[d][100])/2.0;
    }

    write_output(output_filename, &output);


}

void run_benchmark(DerivFunc derivs, double beta_scalar, ScenarioParameters *params)
{
    Matrix<double, NEA, NEA, RowMajor> china_home = Map<Matrix<double, NEA, NEA, RowMajor>>(&China_home[0][0]);
    Matrix<double, NEA, NEA, RowMajor> china_work = Map<Matrix<double, NEA, NEA, RowMajor>>(&China_work[0][0]);
    Matrix<double, NEA, NEA, RowMajor> china_school = Map<Matrix<double, NEA, NEA, RowMajor>>(&China_school[0][0]);
    Matrix<double, NEA, NEA, RowMajor> china_other = Map<Matrix<double, NEA, NEA, RowMajor>>(&China_other[0][0]);

    Matrix<double, NEA, NEA, RowMajor> china_C = china_home + china_work + china_school + china_other;
    Matrix<double, NEA, NEA, RowMajor> china_P;
    for (int i = 0; i < NEA; i++)
    {
        for (int j = 0; j < NEA; j++)
        {
            if ( i == j)
            {
                china_P(i, j) = 1;
            }
            else
            {
                china_P(i, j) = p_vec[i] / p_vec[j];
            }
        }
    }
    double min_r0 = 1.6;
    double max_r0 = 3.0;

    // di = 3, I0 = 200
    double f0 = 0.0002;
    Matrix<double, NEA, NEA, RowMajor> china_P_esc = china_P * f0;;
    Matrix<double, NEA, NEA, RowMajor> china_Cnew = 0.5 *(china_C + china_C.transpose()*china_P_esc) * china_P_esc;
    Matrix<double, NEA, NEA, RowMajor> china_Cesc = china_Cnew/(china_Cnew.eigenvalues().real().maxCoeff());
    run_trial(derivs, params, 3, 200, min_r0, max_r0, china_Cesc, "output/result_b1.csv");

    // di = 3, I0 = 2000
    f0 = 0.002;
    china_P_esc = china_P * f0;;
    china_Cnew = 0.5 *(china_C + china_C.transpose()*china_P_esc) * china_P_esc;
    china_Cesc = china_Cnew/(china_Cnew.eigenvalues().real().maxCoeff());
    run_trial(derivs, params, 3, 2000, min_r0, max_r0, china_Cesc, "output/result_b2.csv");

    // di = 7, I0 = 200
    f0 = 0.0002;
    china_P_esc = china_P * f0;;
    china_Cnew = 0.5 *(china_C + china_C.transpose()*china_P_esc) * china_P_esc;
    china_Cesc = china_Cnew/(china_Cnew.eigenvalues().real().maxCoeff());
    run_trial(derivs, params, 7, 200, min_r0, max_r0, china_Cesc, "output/result_b3.csv");

    // di = 7, I0 = 2000
    f0 = 0.002;
    china_P_esc = china_P * f0;;
    china_Cnew = 0.5 *(china_C + china_C.transpose()*china_P_esc) * china_P_esc;
    china_Cesc = china_Cnew/(china_Cnew.eigenvalues().real().maxCoeff());
    run_trial(derivs, params, 7, 2000, min_r0, max_r0, china_Cesc, "output/result_b4.csv");

} */

int write_output(const char* filename, ScenarioOutput *o)
{
	FILE *result;
	if ((result = fopen(filename, "w")) == NULL) {
		perror("Abertura de arquivo de saida falhou. \n");
		return 1;
	}
    int day = 0;
    double (*YOUT)[NA][MAX_DAYS] = o->YOUT;
    double (*Y_sum)[MAX_DAYS] = o->Y_sum;

	// Print header
	fprintf(result, "S,E,I,R,N,A,C,D,H,L,Qi,Qa,");
	const char* paramNames[] = { "S", "E", "I", "R", "N", "A", "C", "D", "H", "L", "Qi", "Qa" };
	for (int k = 0; k < NEA; k++)
	{
		for (int j = 0; j < (NA - 1); j++)
		{
			fprintf(result, "%s_%d,", paramNames[j], k);
		}
	}
	// Header do Ri mais � direita da tabela
	fprintf(result, "Ri,");
	for (int k = 0; k < NEA; k++)
	{
		fprintf(result, "Ri_%d,", k);
	}
	fprintf(result, "\n");

	// Printa resultado estratificado (dia 0)
    fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,",
              Y_sum[Var::S][day], Y_sum[Var::E][day], Y_sum[Var::I][day], Y_sum[Var::R][day], Y_sum[Var::N][day]
            , Y_sum[Var::A][day], Y_sum[Var::C][day], Y_sum[Var::D][day], Y_sum[Var::H][day], Y_sum[Var::L][day], Y_sum[Var::Qi][day], Y_sum[Var::Qa][day]);
    for (int k = 0; k < NEA; k++)
    {
        fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", YOUT[k][Var::S][day], 
                YOUT[k][Var::E][day], YOUT[k][Var::I][day],
                YOUT[k][Var::R][day], YOUT[k][Var::N][day], 
                YOUT[k][Var::A][day], YOUT[k][Var::C][day], 
			    YOUT[k][Var::D][day], YOUT[k][Var::H][day], YOUT[k][Var::L][day], YOUT[k][Var::Qi][day], YOUT[k][Var::Qa][day]);
    }
	// Printa resultado estratificado do Ri (dia 0)
	fprintf(result, "%lf,", Y_sum[Var::Ri][day]);
	for (int k = 0; k < NEA; k++)
	{
		fprintf(result, "%lf,", YOUT[k][Var::Ri][day]);
	}
    fprintf(result, "\n");

	// Printa resultado estratificado
    for (int i = 1; i < MAX_DAYS; i++)
    {
        day = i;
        fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", 
                  Y_sum[Var::S][day], Y_sum[Var::E][day], Y_sum[Var::I][day], Y_sum[Var::R][day], Y_sum[Var::N][day]
                , Y_sum[Var::A][day], Y_sum[Var::C][day], Y_sum[Var::D][day], Y_sum[Var::H][day], Y_sum[Var::L][day], Y_sum[Var::Qi][day], Y_sum[Var::Qa][day]);
        for (int k = 0; k < NEA; k++)
        {
            fprintf(result, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", YOUT[k][Var::S][day], 
                    YOUT[k][Var::E][day], YOUT[k][Var::I][day],
                    YOUT[k][Var::R][day], YOUT[k][Var::N][day], 
                    YOUT[k][Var::A][day], YOUT[k][Var::C][day], 
				    YOUT[k][Var::D][day], YOUT[k][Var::H][day], YOUT[k][Var::L][day], YOUT[k][Var::Qi][day], YOUT[k][Var::Qa][day]);
        }
		// Printa resultado estratificado do Ri
		fprintf(result, "%lf,", Y_sum[Var::Ri][day]);
		for (int k = 0; k < NEA; k++)
		{
			fprintf(result, "%lf,", YOUT[k][Var::Ri][day]);
		}
        fprintf(result, "\n");
    }
    fclose(result);
	return 0;
}



