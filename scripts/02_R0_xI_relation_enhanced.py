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

# Definitions
I0 = 1
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [0, 6, 8, 9]

y = []
xi_values = np.linspace(0, 0.98, 100)

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

    # First run to generate betas
    y = []
    for xi in xi_values:
        xI_array = np.full(age_strata, xi, dtype=np.float64)
        xA_array = np.zeros(age_strata, dtype=np.float64)
        R0 = R0FromxIxA(scenario, xI_array, xA_array, verbose = False)
        y.append([xi, R0])

    np.savetxt("../../output/results/2_R0_xI_relation/cenario" + scenario + "/itv=" + str(itv) + ".csv", y,
                   delimiter=",", header="x, r_0")
