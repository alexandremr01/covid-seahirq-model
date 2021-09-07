from utils import OutputData, write_testing_parameters, get_day_to_start_testing
from models import SeahirModel
import numpy as np
from model.nextgen import calculateR0, R0FromBetaGama, R0FromxIxA
import os
import itertools
from model.scen_gen import scen_gen, TESTING_PARAMETERS_DIRECT
import pandas as pd
from model.parameters import Parameters
from models import SeahirModel, EpidemicPhases
import matplotlib.pyplot as plt

# Definitions
I0 = 5
scenario='BR'
age_strata = 16
interventions = [0, 6, 8, 9]
default_figsize = [6.4, 4.8]

P_I_Test = 0.5
Pinfec = 0.8

xI = 0.3
xA = 0.3

xI_array = np.full(age_strata, xI, dtype=np.float64)
xA_array = np.full(age_strata, xA, dtype=np.float64)
days = list(range(400))
y = []
fig, axs = plt.subplots(1, 4, figsize=(4*default_figsize[0], default_figsize[1]))

for k, itv in enumerate(interventions):
    print(scenario, ' ', itv)
    parameters = Parameters(age_strata, num_phases=2, testing_parameters_mode=TESTING_PARAMETERS_DIRECT, changeDir=True, scenario=scenario)
    model = SeahirModel(scenario)

    xI_array = np.zeros(age_strata, dtype=np.float64)
    xA_array = np.zeros(age_strata, dtype=np.float64)

    # First run to generate betas
    phases = EpidemicPhases(g0=1, itv0=0)
    phases.add_phase(g=1, start_day=10, itv_level=itv)
    model.run(I0=I0, R0=None, phases=phases, verbose=False, attack=1.0, xI=xI_array, xA=xA_array, testing_parameters = TESTING_PARAMETERS_DIRECT)

    outData = OutputData(scenario, parameters)

    P_Test_max_sint = np.divide(outData.I, outData.N) * xI / P_I_Test
    P_Test_Rand = np.full(400, xI)
    P_Test_Max_infec = np.divide(outData.I*xI + outData.A*xA, outData.N*Pinfec)

    N_Tests_max_sint = P_Test_max_sint * 1e7 # all populations are normalized
    N_Tests_rand = P_Test_Rand * 1e7 # all populations are normalized
    N_Tests_max_infec = P_Test_Max_infec * 1e7 # all populations are normalized

    print("Maximizing sintomatic:", np.sum(N_Tests_max_sint)/1e6)

    print("Maximizing infected:", np.sum(N_Tests_max_infec)/1e6)

    print("Random:", np.sum(N_Tests_rand)/1e6)
    # y = np.concatenate((days, N_Tests_max_sint, N_Tests_rand, N_Tests_max_sint))

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
    axs[k].set_title(name+itvname, {'fontsize': 20})

    # axs[k].set_title(scenario + " - itv=" + str(itv))

    axs[k].set_xlabel("Dia")
    axs[k].plot(days, N_Tests_max_sint, label='Maximizar sintomáticos')
    axs[k].plot(days, N_Tests_max_infec, label='Maximizar infectados')
    axs[k].legend(fontsize=11)

    axs[k].set_xlabel('Dia', fontsize=16)
    axs[k].set_ylabel('Testes', fontsize=16)
    axs[k].tick_params(axis='x', labelsize=14)
    axs[k].tick_params(axis='y', labelsize=14)

plt.tight_layout()
plt.savefig("../../output/results/3_PTest_per_xI/cenario" + scenario + "/result.png")
plt.close()

    # np.savetxt("../../output/results/3_PTest_per_xI/cenario" + scenario + "/itv=" + str(itv) + ".csv", y,
    #                delimiter=",", header="x, r_0")

