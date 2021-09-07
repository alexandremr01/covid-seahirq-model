#pragma once

#include "macros.h"


void rk4(double y[], double dydx[], double x, double h, double yout[],
	void (*derivs)(double, double [], double []));

void rk42D(double y[][NA], double dydx[][NA], double x, double h, double yout[][NA],
	DerivFunc derivs, ScenarioParameters *params);
