from utils import Parameters, OutputData, write_testing_parameters, get_day_to_start_testing, second_derivatives, first_derivatives, doubling_time
from models import SeahirModel
import numpy as np
from plot_test_simulations import plot_fit
import os
import csv

# Código de teste para um xI específico

# Definitions
I0 = 1
uf = 'SP'
parameters = Parameters(uf)
model = SeahirModel(uf)

xA = 0
xI = 0
write_testing_parameters(uf, xI, xA)

# Modelo com intervenção
d0=parameters.days[0]
model.run(I0=I0, R0=parameters.r_0, g=parameters.g, days=parameters.days, testing_parameters=True, itv=[10, 9, 8], verbose=False)
#######################

# Modelo sem intervenção
parameters.reset()
#d0 = get_day_to_start_testing(uf, 10000, parameters)
# d0 = 56
# model.run_ifrs(I0=I0, R0=2.0, zeta=0, tau=parameters.tau, g=[1], days=[d0],
#                    n=2, itv=[0], ifrs=parameters.ifrs, testing_parameters=True, verbose=True)
########################

outData = OutputData(uf, parameters)
maxH = np.max(outData.H)
maxI = np.max(outData.I)
deaths = outData.CD[-1]
#duration = np.max(np.argwhere(outData.I > 1))
print("Segunda derivada por dia: ")
derivatives = first_derivatives(outData.I+outData.Qi)
print(derivatives)
d_virada = np.min(np.argwhere(derivatives < 0.0))
print("Dia 0: ", d0)
print("Dia da virada da derivada: ",d_virada)
print("Half time a partir do dia da virada: ", doubling_time(outData.I+outData.Qi, d_virada, 60))
print("No ultimo dia: "+str(outData.I[int(parameters.days[0])+60]))
print("MaxH: "+str(maxH))
print("MaxI: "+str(maxI))
print("Mortes: "+str(deaths))
#print("Duração: "+str(duration))

h_time = doubling_time(outData.I+outData.Qi, d_virada, 60)
text = '\n'.join(
    [r'Dia 0: %d'%(d0), 'Dia de virada da 1a derivada: %d'%(d_virada), 'Half time a partir da virada: %.2f'%(h_time )])

plot_fit(uf, parameters, d0, xI, text)
