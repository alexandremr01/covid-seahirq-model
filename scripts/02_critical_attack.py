from utils import OutputData, write_testing_parameters, get_day_to_start_testing
import numpy as np
from model import parameters
import cmodels
import os

# Definitions
I0 = 1
scenarios = ['US']
age_strata = 16
interventions = [6, 8, 9]
epsilon = 1e-1
xi_values = np.linspace(0, 0.98, 50)
xa_values = np.linspace(0, 0.98, 50)
xI_array = np.full(age_strata, 0, dtype=np.float64)
xA_array = np.full(age_strata, 0, dtype=np.float64)
attack_values = np.linspace(0.4, 2.5, 40)
    
for scenario in scenarios:
    output_folder = "../output/results/02_critical_attack/scenario%s" % scenario
    os.makedirs(output_folder, exist_ok=True)
    for itv in interventions:
        print(scenario, ' ', itv)
        y = []
        
        phases = parameters.EpidemicPhases(g0=1, itv0=0, xI0=xI_array, xA0=xA_array)
        p = parameters.Parameters(
            scenario=scenario, phases=phases, age_strata=age_strata, model=3, fatality=1.0)

        for xi in xi_values:
            xI_array = np.full(age_strata, xi, dtype=np.float64)
            for xa in xa_values:
                xA_array = np.full(age_strata, xa, dtype=np.float64)
                for attack in attack_values:
                    print('xI=%f xA=%f attack=%f' % (xi, xa, attack))
                    phases = parameters.EpidemicPhases(g0=1, itv0=itv, xI0=xI_array, xA0=xA_array)
                    R0 = p.update_phases(phases, attack=attack)
                    print(R0)
                    if abs(R0 - 1.0) < epsilon:
                        y.append([xi, xa, attack])
                        break

        output_path = "%s/itv=%s.csv" % (output_folder, str(itv))
        np.savetxt(output_path, y, delimiter=",", header="xI, xA, attack")
