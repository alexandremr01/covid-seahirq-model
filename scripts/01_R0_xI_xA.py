from utils import OutputData, write_testing_parameters, get_day_to_start_testing
import numpy as np
from model import parameters
import cmodels
import os

# Definitions
I0 = 1
scenarios = ['BR', 'NG', 'US', 'GR']
age_strata = 16
interventions = [8, 9]
xi_values = np.linspace(0, 0.98, 100)
xa_values = np.linspace(0, 0.98, 100)
xI_array = np.full(age_strata, 0, dtype=np.float64)
xA_array = np.full(age_strata, 0, dtype=np.float64)
    
for scenario in scenarios:
    output_folder = "../output/results/01_R0_xI_xa/scenario%s" % scenario
    os.makedirs(output_folder, exist_ok=True)
    for itv in interventions:
        print(scenario, ' ', itv)
        y = []
        
        phases = parameters.EpidemicPhases(g0=1, itv0=0, xI0=xI_array, xA0=xA_array)
        p = parameters.Parameters(
            scenario=scenario, phases=phases, age_strata=age_strata, model=3, fatality=1.0)

        for xi in xi_values:
            print (xi)
            xI_array = np.full(age_strata, xi, dtype=np.float64)
            for xa in xa_values:
                xA_array = np.full(age_strata, xa, dtype=np.float64)
                phases = parameters.EpidemicPhases(g0=1, itv0=itv, xI0=xI_array, xA0=xA_array)
                R0 = p.update_phases(phases)
                y.append([xi, xa, R0])

        output_path = "%s/itv=%s.csv" % (output_folder, str(itv))
        np.savetxt(output_path, y, delimiter=",", header="xI, xA, R0")
