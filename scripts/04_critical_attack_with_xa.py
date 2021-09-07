from utils import OutputData, write_testing_parameters, get_day_to_start_testing
from models import SeahirModel
import numpy as np
from model.nextgen import calculateR0, R0FromBetaGama, R0FromxIxA, NextGenScenario
import os
import itertools
from model.scen_gen import scen_gen, TESTING_PARAMETERS_DIRECT
import pandas as pd
from model.parameters import Parameters
from models import SeahirModel, EpidemicPhases

# Definitions
I0 = 1
epsilon = 1e-1
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [0, 6, 8, 9]

y = []
xi_values = np.linspace(0, 0.98, 100)
attack_values = np.linspace(0, 2.0, 100)

for scenario, itv in list(itertools.product(scenarios, interventions)):
    print(scenario, ' ', itv)
    parameters = Parameters(age_strata, num_phases=2, testing_parameters_mode=TESTING_PARAMETERS_DIRECT, changeDir=True, scenario=scenario)
    model = SeahirModel(scenario)

    xI_array = np.zeros(age_strata, dtype=np.float64)
    xA_array = np.zeros(age_strata, dtype=np.float64)

    # First run to generate betas
    phases = EpidemicPhases(g0=1, itv0=0)
    phases.add_phase(g=1, start_day=10, itv_level=itv)
    model.run(I0=I0, R0=None, phases=phases, verbose=False, attack=1.0, xI=xI_array, xA=xA_array, testing_parameters = TESTING_PARAMETERS_DIRECT)
    nextGenScenario = NextGenScenario(scenario, changeDir=True)

    # First run to generate betas
    y = []
    for xi in xi_values:
        print("xI=",xi)
        for xa in xi_values:
            xI_array = np.full(age_strata, xi, dtype=np.float64)
            xA_array = np.full(age_strata, xa, dtype=np.float64)
            for attack in attack_values:
                R0 = R0FromxIxA(scenario, xI_array, xA_array, changeDir=True, attack=attack, verbose = False, nextGenScenario=nextGenScenario)
                if abs(R0 - 1.0) < epsilon:
                        y.append([xi, xa, attack])
                        break

    np.savetxt("../output/results/4_critical_attack/cenario" + scenario + "/xa_itv=" + str(itv) + ".csv", y,
                   delimiter=",", header="x, attack")
