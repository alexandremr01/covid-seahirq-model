from utils import OutputData, write_testing_parameters, get_day_to_start_testing, second_derivatives, first_derivatives, doubling_time
from models import SeahirModel, EpidemicPhases
import numpy as np
from plot_test_simulations import plot_fit
import os
import csv
from model.parameters import Parameters
from model.scen_gen import scen_gen, TESTING_PARAMETERS_FROM_EPIDEMIOLOGY


# Código de teste para um xI específico

# Definitions
I0 = 1
uf = 'GR'
age_strata=16
parameters = Parameters(age_strata, num_phases=2, testing_parameters_mode=TESTING_PARAMETERS_FROM_EPIDEMIOLOGY, changeDir=True, scenario=uf)
model = SeahirModel(uf)

xA = 0
xI = 0
write_testing_parameters(uf, xI, xA)

# Modelo com intervenção
model.run(I0=I0, R0=2.0, phases=EpidemicPhases(g0=1, itv0=0), 
                  verbose=False, attack=None, testing_parameters = 1)
#######################
d0=0

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
print("Infectados no ultimo dia: "+str(outData.I[-1]))

print("MaxH: "+str(maxH))
print("MaxI: "+str(maxI))
print("Mortes: "+str(deaths))
#print("Duração: "+str(duration))

h_time = doubling_time(outData.I+outData.Qi, d_virada, 60)
text = '\n'.join(
    [r'Dia 0: %d'%(d0), 'Dia de virada da 1a derivada: %d'%(d_virada), 'Half time a partir da virada: %.2f'%(h_time )])

plot_fit(uf, parameters, d0, xI, text)
