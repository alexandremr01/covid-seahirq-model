import numpy as np
import matplotlib.pyplot as plt
from utils import second_derivatives, Parameters, OutputData
from models import SeahirModel
import itertools

I0 = 1
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [[0], [6], [8], [9]]

for scenario, itv in list(itertools.product(scenarios, interventions)):
    y = np.loadtxt("../output/results/4_critical_attack/cenario" + scenario + "/itv=" + str(itv) + ".csv", delimiter=',')
    attack = y[:, 1]
    xI = y[:, 0]
    plt.scatter(xI, attack)
    plt.xlabel('xI')
    plt.ylabel('attack')
    plt.title('Susceptibility por xI com itv='+str(itv) + " for "+ scenario)
    plt.savefig("../output/results/4_critical_attack/cenario" + scenario + "/itv=" + str(itv) + ".png")
    #plt.show()
    plt.close()