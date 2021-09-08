import cmodels
import numpy as np
from modelv2 import parameters

model = 3
age_strata=16

phases = parameters.EpidemicPhases(g0=1, itv0=0, xI0=np.full(age_strata, 0, dtype=np.float64), xA0=np.full(age_strata, 0, dtype=np.float64))

p = parameters.Parameters(scenario='BR', phases=phases, age_strata=age_strata, model=model, fatality=1.0)

out = cmodels.model(model, p.initial, p.cparameters, p.dynamic_parameters)

# print(p.initial)
# print(out.Y_sum)