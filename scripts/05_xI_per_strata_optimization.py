from utils import Parameters, OutputData, write_testing_parameters, doubling_time, age_strata, first_derivatives, get_day_to_start_testing
from models import SeahirModel
import numpy as np
from plot_test_simulations import plot_fit
import os
import csv
from scipy import optimize


# Definitions
I0 = 1
age_strata = 16
uf = 'SP'
I_to_start_testing = 10
parameters = Parameters(uf)
model = SeahirModel(uf)
xA_array = np.zeros(age_strata, dtype=np.float64)

# Simulações para testes com primeira intervenção
def cost(xI):
    model.run(I0=I0, R0=2.0, zeta=0.01, tau=2.9, g=[1], days=[d_0], n=2,
                      itv=[0], xI=xI, xA=xA_array, verbose=False)
    outData = OutputData(uf, parameters)

    H0 = 14
    derivatives = first_derivatives(outData.I + outData.Qi)
    # print(derivatives)
    # input()
    derivatives_aux = derivatives[d_0:]
    d_max = np.min(np.argwhere(derivatives_aux < 0.01))            # primeiro dia de derivada 0 apos o inicio
    half_life = doubling_time(outData.I + outData.Qi, d_max+d_0, 60)
    # d_max = d_virada - d_0 # tempo até virada
    cost_xI = np.sum(np.power(xI, 2))/16

    J = H0*cost_xI + half_life*half_life + d_max*d_max
    print("Cost:", J, "; cost_xI=", cost_xI, " half_life=",half_life, " d=",d_max)
    return J

model.run(I0=I0, R0=2.0, zeta=0.01, tau=2.9, g=[1], days=[0], n=2,
                      itv=[0], verbose=False)
d_0 = get_day_to_start_testing(uf, I_to_start_testing, parameters)
bounds = [[0.0, 1.0]] * 16 #mesmos bounds para os 16 parâmetros
n_iteration_ref = 20
pop_size = 30
ret_ref = optimize.differential_evolution(cost, bounds,
        args=(), disp=True, maxiter=n_iteration_ref,
        popsize=pop_size, seed=1234)