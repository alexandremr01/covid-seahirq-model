import cmodels
import numpy as np
from model import parameters
import os

def death_day_relation():
    I0 = 1
    scenario = 'BR'
    age_strata = 16
    xi_set = [ 0.9]

    nullXI = np.full(age_strata, 0, dtype=np.float64)

    phases = parameters.EpidemicPhases(
        g0=1, itv0=0, xI0=nullXI, xA0=nullXI)
    p = parameters.Parameters(
        scenario='BR', phases=phases, age_strata=age_strata, model=3, fatality=1.0)
    output_folder = "../output/results/03_death_day_relation/scenario%s" % scenario
    os.makedirs(output_folder, exist_ok=True)

    for xI in xi_set:
        for xA in xi_set:
            print(scenario, ' ', xI)

            y = []
            xI_array = np.full(age_strata, xI, dtype=np.float64)
            xA_array = np.full(age_strata, xA, dtype=np.float64)

            for d in range(400):
                print("Running day ", d)
                # Runs start testing on day d
                phases = parameters.EpidemicPhases(
                    g0=1, itv0=0, xI0=nullXI, xA0=nullXI)
                phases.add_phase(g=1, start_day=d, itv_level=0,
                                 xI=xI_array, xA=xA_array)
                p.update_phases(phases)
                print(p.initial)
                out = cmodels.model(
                    3, p.initial, p.cparameters, p.dynamic_parameters, False)
                # print(out.Y_sum[-1])
                # print(out.YOUT)
                deaths = out.Y_sum[-1][10]
                print(out.Y_sum[-1][10])
                y.append([d, deaths])

            output_path = output_folder + "/xI=%.2f_xA=%.2f.csv" % (xI, xA)
            np.savetxt(output_path, y, delimiter=",", header="day, deaths")


if __name__ == '__main__':
    death_day_relation()
