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
scenario = 'BR'
age_strata = 16
xi_set = [0.2, 0.5, 0.7, 0.9]

model = SeahirModel(scenario)
for xI in xi_set:
    for xA in xi_set:
        print(scenario, ' ', xI)

        x = list(range(400))
        y = []
        xI_array = np.full(age_strata, xI, dtype=np.float64)
        xA_array = np.full(age_strata, xA, dtype=np.float64)

        for d in range(400):
            print("Running day ", d)
            # Runs start testing on day d
            phases = EpidemicPhases(g0=1, itv0=0)
            phases.add_phase(g=1, start_day=d, itv_level=0)
            model.run(I0=I0, R0=None, phases=phases, verbose=False, attack=1.0, xI=xI_array, xA=xA_array, testing_parameters = TESTING_PARAMETERS_DIRECT)

            # Store number of deaths
            outData = OutputData(scenario)
            deaths = outData.CD[-1]
            print(deaths)
            y.append([d, deaths])

        np.savetxt("../output/results/1_test_start/cenario" + scenario + "/data_xI=" + str(xI) + "_xA="+ str(xA) + ".csv", y,
                delimiter=",", header="day, deaths")
