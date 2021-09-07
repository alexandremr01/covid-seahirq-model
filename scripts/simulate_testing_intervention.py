from utils import Parameters, OutputData, write_testing_parameters
from models import SeahirModel
import numpy as np
from plot_test_simulations import plot_fit
import os
import csv

# Simulações para testes com primeira intervenção

# Definitions
I0 = 1
uf = 'SP'
I_to_start_testing = 10
parameters = Parameters(uf)
model = SeahirModel(uf)

I_set = [10, 100, 1000, 10000]
#Tabela para guardar, para cada I: maxH, mortos, duracao por xI
result_tables = np.zeros([4, 100, 4], dtype=np.float64)
# First run - discover day with X infected
# parameters.reset()
model.run_ifrs(I0=I0, R0=parameters.r_0, zeta=parameters.zeta, tau=parameters.tau, g=parameters.g, days=parameters.days,
               itv=[10, 9, 8], ifrs=parameters.ifrs, verbose=False)
succesfull_xI = 0
xA = 0

for j, xI in enumerate(np.arange(0, 1.0, 0.01)):
    print("Analyzing " + str(xI))
    write_testing_parameters(uf, xI, xA)
    model.run_ifrs(I0=I0, R0=parameters.r_0, zeta=parameters.zeta, tau=parameters.tau, g=parameters.g,
               days=parameters.days, testing_parameters=True, itv=[10, 9, 8], ifrs=parameters.ifrs, verbose=False)
    outData = OutputData(uf, parameters)
    maxH = np.max(outData.H)
    deaths = outData.CD[-1]
    duration = np.max(np.argwhere(outData.I > 1))
    print("No ultimo dia: "+str(outData.I[int(parameters.days[0])+60]))
    if outData.I[int(parameters.days[0])+30] < 1:                          # Em 30 dias, consegue parar a doença?
        succesfull_xI = xI
        break
print("xI minimo para fim da pandemia: " + str(succesfull_xI))
plot_fit(uf, parameters, 0)
