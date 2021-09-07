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
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [0, 6, 8, 9]

P_I_Test = 0.1
target_xI = 0.5

xI_array = np.full(age_strata, target_xI, dtype=np.float64)
xA_array = np.zeros(age_strata, dtype=np.float64)
days = list(range(400))
y = []
for scenario, itv in list(itertools.product(scenarios, interventions)):
    print(scenario, ' ', itv)
    parameters = Parameters(age_strata, num_phases=2, testing_parameters_mode=TESTING_PARAMETERS_DIRECT, changeDir=True, scenario=scenario)
    model = SeahirModel(scenario)

    xI_array = np.zeros(age_strata, dtype=np.float64)
    xA_array = np.zeros(age_strata, dtype=np.float64)

    # First run to generate betas
    phases = EpidemicPhases(g0=1, itv0=0)
    phases.add_phase(g=1, start_day=50, itv_level=itv)
    model.run(I0=I0, R0=None, phases=phases, verbose=False, attack=1.0, xI=xI_array, xA=xA_array, testing_parameters = TESTING_PARAMETERS_DIRECT)

    outData = OutputData(scenario, parameters)
    P_Test = np.divide(outData.I, outData.N) * target_xI / P_I_Test

    N_Tests = P_Test * 1e7 # all populations are normalized
    y = np.concatenate((days, N_Tests))

    plt.title(scenario + " - itv=" + str(itv))
    plt.xlabel("Dia")
    plt.ylabel("# testes")
    plt.scatter(days, N_Tests)
    plt.savefig("../output/results/3_PTest_per_xI/cenario" + scenario + "/itv=" + str(itv) + ".png")
    plt.close()

    np.savetxt("../output/results/3_PTest_per_xI/cenario" + scenario + "/itv=" + str(itv) + ".csv", y,
                   delimiter=",", header="x, r_0")

