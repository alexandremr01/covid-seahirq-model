from utils import OutputData, write_testing_parameters, get_day_to_start_testing
import numpy as np
import os
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from model import parameters
import cmodels 
# Definitions
I0 = 5
scenarios=['BR', 'GR', 'US', 'NG']
age_strata = 16
interventions = [0, 6, 8, 9]
default_figsize = [6.4, 4.8]

P_I_Test = 0.5
Pinfec = 0.5

xI = 0.3
xA = 0.3

xI_array = np.full(age_strata, xI, dtype=np.float64)
xA_array = np.full(age_strata, xA, dtype=np.float64)
days = list(range(400))
y = []
# fig, axs = plt.subplots(1, 4, figsize=(4*default_figsize[0], default_figsize[1]))

default_figsize = [6.4, 4.8]
fig, axs = plt.subplots(4, 4, figsize=(4*default_figsize[0], 4*default_figsize[1]))

for m, scenario in enumerate(scenarios):
    for k, itv in enumerate(interventions):
        print(scenario, ' ', itv)
        phases = parameters.EpidemicPhases(g0=1, itv0=itv, xI0=xI_array, xA0=xA_array)
        p = parameters.Parameters(scenario=scenario, phases=phases, age_strata=age_strata, model=3, fatality=1.0)
        out = cmodels.model(3, p.initial, p.cparameters, p.dynamic_parameters, False)

        I = np.array([out.Y_sum[x][4] for x in range(len(out.Y_sum))])
        A = np.array([out.Y_sum[x][2] for x in range(len(out.Y_sum))])
        N = np.array([out.Y_sum[x][9] for x in range(len(out.Y_sum))])

        P_Test_max_sint = np.divide(I, N) * xI / P_I_Test
        P_Test_Rand = np.full(400, xI)
        P_Test_Max_infec = np.divide(I*xI + A*xA, N*Pinfec)

        N_Tests_max_sint = P_Test_max_sint * 1e7 # all populations are normalized
        N_Tests_rand = P_Test_Rand * 1e7 # all populations are normalized
        N_Tests_max_infec = P_Test_Max_infec * 1e7 # all populations are normalized

        print("Maximizing sintomatic: %f" % (np.sum(N_Tests_max_sint)/1e6))
        print("Maximizing infected: %f" % (np.sum(N_Tests_max_infec)/1e6))
        print("Random: %f" % (np.sum(N_Tests_rand)/1e6))

        if scenario == 'BR':
            name = "Brasil"
        elif scenario == 'GR':
            name = "Alemanha"
        elif scenario == 'NG':
            name = "Nigéria"
        elif scenario == 'US':
            name = "EUA"

        if itv == 0:
            itvname = ' - Sem intervenção'
        elif itv == 6:
            itvname = ' - Isolamento Vertical'
        elif itv == 8:
            itvname = ' - Distanciamento social'
        elif itv == 9:
            itvname = ' - Com restrições no trabalho'
        axs[m, k].set_title(name+itvname, {'fontsize': 20})

        # axs[k].set_title(scenario + " - itv=" + str(itv))

        axs[m, k].set_xlabel("Dia")
        axs[m, k].plot(days, N_Tests_max_sint, label='Maximizar sintomáticos')
        axs[m, k].plot(days, N_Tests_max_infec, label='Maximizar infectados')
        axs[m, k].legend(fontsize=11)

        axs[m, k].set_xlabel('Dia', fontsize=16)
        axs[m, k].set_ylabel('Testes', fontsize=16)
        axs[m, k].tick_params(axis='x', labelsize=14)
        axs[m, k].tick_params(axis='y', labelsize=14)

plt.tight_layout()
plt.savefig("../output/results/04_PTest_per_xI/result.png")
plt.close()