#pragma once

#include "macros.h"
#include "scenario.h"



// A soma será a da linha i da matrix f(Dij)
double sum_beta(int k, double y[][NA], ScenarioParameters *p);

// A soma será a da linha i da matrix f(Dij)
double sum_alpha_beta(int k, double y[][NA], ScenarioParameters *p);

double sum_exp_beta(int k, double y[][NA], ScenarioParameters *p);

void derivs_sir(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p);

void derivs_seir(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p);

void derivs_seair(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p);

void derivs_seahirq(double t, double y[][NA], double dydt[][NA], ScenarioParameters *p) ;

void update_params(double y[][NA], ScenarioParameters *p);

DerivFunc get_model(Model::Enum model);


