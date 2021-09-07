from utils import OutputData, write_testing_parameters, get_day_to_start_testing
from models import SeahirModel, EpidemicPhases
import numpy as np
from model.nextgen import calculateR0, R0FromBetaGama, R0FromxIxA
import os
import itertools
from model.scen_gen import scen_gen, TESTING_PARAMETERS_FROM_EPIDEMIOLOGY
import pandas as pd
from model.parameters import Parameters

# Definitions
I0 = 5
scenarios = ['BR']
age_strata = 16
interventions = [[0]]

def calculateForce(outData, parameters):
    df = pd.read_csv('../../input/cenarios/cenarioBR/beta_gama.csv', sep=',')
    betas = df.values[0:age_strata, 17 + 2 * age_strata:17 + 3 * age_strata]
    print(np.sum(outData.n[-1]))
    force = np.dot(outData.y[-1], betas) + np.dot(np.multiply(outData.a[-1], parameters.alpha), betas) + 0.4*np.dot(outData.e[-1], betas)
    force = np.divide(force, outData.n[-1])
    return force

# Para cada beta, gera R0, Iinf e plota
for scenario, itv in list(itertools.product(scenarios, interventions)):
    beta_values = np.linspace(0, 1.5, 100)

    
    parameters = Parameters(age_strata, num_phases=2, testing_parameters_mode=TESTING_PARAMETERS_FROM_EPIDEMIOLOGY, changeDir=True, scenario=scenario)
    model = SeahirModel(scenario)

    y = []
    phases = EpidemicPhases(g0=1, itv0=0)

    for beta in beta_values:
        model.run(I0=I0, R0=beta, phases=phases, verbose=False, attack=1.0, testing_parameters = TESTING_PARAMETERS_FROM_EPIDEMIOLOGY)
        outData = OutputData(scenario, parameters)
        Iinf = outData.I[-1]
        deaths = outData.CD[-1]
        Iacc = np.add.accumulate(outData.I, axis=0)
        force = calculateForce(outData, parameters)
        R0 = R0FromBetaGama(scenario, verbose=True)
        y.append([beta, R0, Iinf, Iacc[-1], deaths, np.average(force)])

    np.savetxt("../../output/results/bifurcation_analysis/cenario" + scenario + "/itv=" + str(itv) + ".csv", y,
                   delimiter=",", header="x, r_0")

