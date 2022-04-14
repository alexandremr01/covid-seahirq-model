#include "models.h"
#include <iostream>

// A soma será a da linha i da matrix f(Dij)
double sum_beta(int k, double y[][NA], ScenarioParameters *p) {
    double population = 0;
    for (int i = 0; i < NEA; ++i)
    {
        population+=  y[i][Var::N];
    }

    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum +=  p->beta[p->day][k][j] * (y[k][Var::S]/population) * y[j][Var::I];
    }

    return sum;
}

// A soma será a da linha i da matrix f(Dij)
double sum_alpha_beta(int k, double y[][NA], ScenarioParameters *p) {
    double population = 0;
    for (int i = 0; i < NEA; ++i)
    {
        population+=  y[i][Var::N];
    }

    double sum = 0.0;

    for (int j = 0; j < NEA; ++j)
    {
        sum +=  p->alpha[j] * p->beta[p->day][k][j] * (y[k][Var::S]/population) * y[j][Var::A];
    }

    return sum;
}

double sum_exp_beta(int k, double y[][NA], ScenarioParameters *p) {
	double population = 0;
	for (int i = 0; i < NEA; ++i)
	{
		population += y[i][Var::N];
	}

	double sum = 0.0;

	for (int j = 0; j < NEA; ++j)
	{
		sum += p->ksi[j] * p->beta[p->day][k][j] * (y[k][Var::S] / population) * y[j][Var::E];
	}

	return sum;
}


void derivs_sir(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p)
{
    for (int k = 0; k < NEA; ++k)
    {
        double beta_row_sum = sum_beta(k, y, p);
        dydt[k][Var::S] = p->Lambda[k]*y[k][Var::N] - p->mu_eq[k] * y[k][Var::S] - beta_row_sum;
        dydt[k][Var::I] = beta_row_sum - (p->mu_eq[k] + p->mu_cov[k] + p->gama[p->day][k]) * y[k][Var::I];
        dydt[k][Var::R] = p->gama[p->day][k] * y[k][Var::I] - p->mu_eq[k] * y[k][Var::R];
        dydt[k][Var::N] = p->Lambda[k]*y[k][Var::N] - p->mu_eq[k] * y[k][Var::N] - p->mu_cov[k] * y[k][Var::I];
        dydt[k][Var::C] = p->mu_cov[k] * y[k][Var::I];
        dydt[k][Var::H] = p->phi[k] * beta_row_sum - p->gama[p->day][k] * y[k][Var::H];
		dydt[k][Var::L] = p->eta[k] * beta_row_sum - p->gama[p->day][k] * (1 + p->Tlc[k] / (1 - p->Tlc[k])) * y[k][Var::L];

		dydt[k][Var::E] = 0;
		dydt[k][Var::A] = 0;
		dydt[k][Var::Qi] = 0;
		dydt[k][Var::Qa] = 0;
		dydt[k][Var::Ri] = 0;
		dydt[k][Var::D] = 0;
    }
}

void derivs_seir(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p)
{
    for (int k = 0; k < NEA; ++k)
    {
        double beta_row_sum = sum_beta(k, y, p);
        dydt[k][Var::S] = p->Lambda[k]*y[k][Var::N] - p->mu_eq[k] * y[k][Var::S] - beta_row_sum;
        dydt[k][Var::E] = beta_row_sum - (p->mu_eq[k] + p->a[k]) * y[k][Var::E];
        dydt[k][Var::I] = p->a[k] * y[k][Var::E] - (p->mu_eq[k] + p->mu_cov[k] + p->gama[p->day][k]) * y[k][Var::I];
        dydt[k][Var::R] = p->gama[p->day][k] * y[k][Var::I] - p->mu_eq[k] * y[k][Var::R];
        dydt[k][Var::N] = p->Lambda[k]*y[k][Var::N] - p->mu_eq[k] * y[k][Var::N] - p->mu_cov[k] * y[k][Var::I];
        dydt[k][Var::C] = p->mu_cov[k] * y[k][Var::I];
        dydt[k][Var::H] = p->phi[k] * beta_row_sum - p->gama[p->day][k] * y[k][Var::H];
		dydt[k][Var::L] = p->eta[k] * beta_row_sum - p->gama[p->day][k] * (1 + p->Tlc[k] / (1 - p->Tlc[k])) * y[k][Var::L];

		dydt[k][Var::A] = 0;
		dydt[k][Var::Qi] = 0;
		dydt[k][Var::Qa] = 0;
		dydt[k][Var::Ri] = 0;
		dydt[k][Var::D] = 0;
    }
}

void derivs_seair(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p)
{
    for (int k = 0; k < NEA; ++k)
    {
        double beta_row_sum = sum_beta(k, y, p);
        double alpha_beta_row_sum = sum_alpha_beta(k, y, p);
        dydt[k][Var::S] = p->Lambda[k]*y[k][Var::N] - p->mu_eq[k] * y[k][Var::S] - alpha_beta_row_sum - beta_row_sum;
        dydt[k][Var::E] = beta_row_sum + alpha_beta_row_sum - (p->mu_eq[k] + p->a[k]) * y[k][Var::E];        
        dydt[k][Var::I] = p->a[k] * p->rho[k] * y[k][Var::E] - (p->mu_eq[k] + p->mu_cov[k] + p->gama[p->day][k]) * y[k][Var::I];
		dydt[k][Var::A] = p->a[k] * (1.0 - p->rho[k]) * y[k][Var::E] - (p->mu_eq[k] + p->gama_RA[k]) * y[k][Var::A];
		dydt[k][Var::R] = p->gama[p->day][k] * y[k][Var::I] + p->gama_RA[k] * y[k][Var::A] - p->mu_eq[k] * y[k][Var::R];
		dydt[k][Var::Ri] = p->gama[p->day][k] * y[k][Var::I];
        dydt[k][Var::N] = p->Lambda[k]*y[k][Var::N] - p->mu_eq[k] * y[k][Var::N] - p->mu_cov[k] * y[k][Var::I];
        dydt[k][Var::C] = p->mu_cov[k] * y[k][Var::I];
        dydt[k][Var::H] = p->phi[k] * beta_row_sum - p->gama[p->day][k] * y[k][Var::H];
		dydt[k][Var::L] = p->eta[k] * beta_row_sum - p->gama[p->day][k] * (1 + p->Tlc[k] / (1 - p->Tlc[k])) * y[k][Var::L];


		dydt[k][Var::D] = 0;
		dydt[k][Var::Qi] = 0;
		dydt[k][Var::Qa] = 0;
    }
}

void derivs_seahirq(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p) 
{
  for (int k = 0; k < NEA; ++k)
  {
      double beta_row_sum = sum_beta(k, y, p);
      double alpha_beta_row_sum = sum_alpha_beta(k, y, p);
	  double exp_beta_row_sum = sum_exp_beta(k, y, p);
	  dydt[k][Var::S] = p->Lambda[k]*y[k][Var::N] -  p->mu_eq[k] * y[k][Var::S] - alpha_beta_row_sum - beta_row_sum - exp_beta_row_sum;
      dydt[k][Var::E] = beta_row_sum + alpha_beta_row_sum + exp_beta_row_sum  - (p->mu_eq[k] + p->a[k]) * y[k][Var::E];
	  dydt[k][Var::A] = p->a[k] * (1.0 - p->rho[k]) * y[k][Var::E] - (p->mu_eq[k] + p->gama_RA[k] + p->gama_QA[p->day][k]) * y[k][Var::A];
	  dydt[k][Var::I] = p->a[k] * p->rho[k] * y[k][Var::E] - (p->theta[k] * p->mu_cov[k] + p->mu_eq[k] + p->gama_H[k] + p->gama_RI[k] + p->gama_QI[p->day][k]) * y[k][Var::I];
	  dydt[k][Var::Qi] = p->gama_QI[p->day][k] * y[k][Var::I] - (p->gama_HQI[k] + p->mu_eq[k] + p->gama_RQI[k]) * y[k][Var::Qi];
	  dydt[k][Var::Qa] = p->gama_QA[p->day][k] * y[k][Var::A] - (p->gama_RQA[k] + p->mu_eq[k]) * y[k][Var::Qa];
      dydt[k][Var::H] = p->gama_H[k] * y[k][Var::I] + p->gama_HQI[k] * y[k][Var::Qi] - (p->gama_HR[k] + p->mu_eq[k] + p->mu_cov[k]) * y[k][Var::H];
      dydt[k][Var::R] = p->gama_RI[k] * y[k][Var::I] + p->gama_RA[k] * y[k][Var::A] + p->gama_RQA[k] * y[k][Var::Qa] + p->gama_RQI[k] * y[k][Var::Qi] + p->gama_HR[k] * y[k][Var::H] - p->mu_eq[k] * y[k][Var::R];
	  dydt[k][Var::Ri] = p->gama_RI[k] * y[k][Var::I] + p->gama_RQI[k] * y[k][Var::Qi] + p->gama_HR[k] * y[k][Var::H];
	  dydt[k][Var::C] = p->mu_cov[k] * y[k][Var::H];
	  dydt[k][Var::D] = p->theta[k] * p->mu_cov[k] *  y[k][Var::I];
	  dydt[k][Var::N] = p->Lambda[k] * y[k][Var::N] - p->mu_eq[k] * y[k][Var::N] - p->mu_cov[k] * y[k][Var::H] - p->theta[k] * p->mu_cov[k] * y[k][Var::I];
  }
}

void update_params(double y[][NA], ScenarioParameters *p) {
    if (p->model == Model::SEAHIRQ) 
    {
        for (int k = 0; k < NEA; ++k)
        {
			//p->mu_cov[k] = p->gama_HQI[k] * p->Tlc[k] / (1 - p->Tlc[k]);
			y[k][Var::L] = (p->eta[k] / p->phi[k]) * y[k][Var::H];
		    //y[k][Var::L] = (p->eta[k] / p->phi[k]) * y[k][Var::H];
            //y[k][Var::N] = y[k][Var::S] + y[k][Var::E] + y[k][Var::A] + y[k][Var::H]
            //             + y[k][Var::I] + y[k][Var::R] + y[k][Var::Qi] + y[k][Var::Qa];
        }
    }
    /*else 
    {
        for (int k = 0; k < NEA; ++k)
        {
            if (y[k][Var::C] > 0 && y[k][Var::R] > 0)
            {
                p->mu_cov[k] =  p->gama[p->day][k] * y[k][Var::C] / y[k][Var::R];
            }
        }
    }*/
}


DerivFunc get_model(Model::Enum model) {
    DerivFunc f[] = {&derivs_sir, &derivs_seir, &derivs_seair, &derivs_seahirq};
    return f[model];
}

