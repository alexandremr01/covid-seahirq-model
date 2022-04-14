#include "solvers.h"


void rk4(double y[], double dydx[], double x, double h, double yout[],
	void (*derivs)(double, double [], double []))
{
	double xh,hh,h6;

    double dym[NA];
    double dyt[NA];
    double yt[NA];

	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;

	for (int i = 0; i < NA; i++) 
    {
        yt[i] = y[i] + hh * dydx[i];
    }
	(*derivs)(xh,yt,dyt);

	for (int i = 0; i < NA; i++) 
    {
        yt[i] = y[i] + hh * dyt[i];
    }
	(*derivs)(xh,yt,dym);

	for (int i = 0; i < NA; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);

	for (int i = 0; i < NA; i++)
    {
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    }

}

void rk42D(double y[][NA], double dydx[][NA], double x, double h, double yout[][NA],
	DerivFunc derivs, ScenarioParameters *params)
{
	double xh,hh,h6;

    double dym[NEA][NA] = {{0}};
    double dyt[NEA][NA] = {{0}};
    double yt[NEA][NA] = {{0}};

	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;

    for (int i = 0; i < NEA; ++i) 
    {
        for (int j = 0; j < NA; ++j)
        {
            yt[i][j] = y[i][j] + hh * dydx[i][j];
        }
    }
	(*derivs)(xh,yt,dyt, params);

    for (int i = 0; i < NEA; ++i) 
    {
        for (int j = 0; j < NA; ++j)
        {
            yt[i][j] = y[i][j]+ hh * dyt[i][j];
        }
    }
	(*derivs)(xh,yt,dym, params);

    for (int i = 0; i < NEA; ++i) 
    {
        for (int j = 0; j < NA; ++j)
        {
            yt[i][j] = y[i][j] + h * dym[i][j];
            dym[i][j] += dyt[i][j];
        }
    }
	(*derivs)(x+h,yt,dyt, params);

    for (int i = 0; i < NEA; ++i) 
    {
        for (int j = 0; j < NA; ++j)
        {
            yout[i][j] = y[i][j] + h6 * (dydx[i][j] + dyt[i][j] + 2.0 * dym[i][j]);
        }
    }
}



