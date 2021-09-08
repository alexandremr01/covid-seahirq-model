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

import cmodels
import numpy as np
from modelv2 import parameters
import time

model = 3


# Definitions
I0 = 1
scenario = 'BR'
age_strata = 16
xi_set = [0.2, 0.5, 0.7, 0.9]

model = SeahirModel(scenario)
start = time.time()
prev_time = start
all_times = []

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
            phases = parameters.EpidemicPhases(g0=1, itv0=0, xI0=xI_array, xA0=xA_array)
            phases.add_phase(g=1, start_day=d, itv_level=0, xI=xI_array, xA=xA_array)

            p = parameters.Parameters(scenario='BR', phases=phases, age_strata=age_strata, model=3, fatality=1.0)

            out = cmodels.model(3, p.initial, p.cparameters, p.dynamic_parameters)

            deaths = out.Y_sum[-1][10]
            print(out.Y_sum[-1][10])
            y.append([d, deaths])

            new_time = time.time()
            all_times.append(new_time - prev_time)
            prev_time=new_time


        np.savetxt("../output/results/1_test_start/cenario" + scenario + "/data_xI=" + str(xI) + "_xA="+ str(xA) + ".csv", y,
                delimiter=",", header="day, deaths")
end=time.time()
print("Total time")
print(end-start)
all_np_times = np.array(all_times)
print(np.average(all_np_times))
print(np.std(all_np_times))