import numpy as np
import matplotlib.pyplot as plt
from utils import second_derivatives, OutputData
from models import SeahirModel
import itertools

I0 = 1
scenarios = ['BR']
age_strata = 16
interventions = [[0]]

for scenario, itv in list(itertools.product(scenarios, interventions)):
    y = np.loadtxt("../output/results/bifurcation_analysis/cenario" + scenario + "/itv=" + str(itv) + ".csv", delimiter=',')
    Iinf = y[:, 2]
    Iacc = y[:, 3]
    death = y[:, 4]
    force = y[:, 5]
    R0 = y[:, 1]
    plt.plot(R0, force)
    plt.xlabel('R0')
    plt.ylabel('Force')
    plt.title('Gráfico de Bifurcação')
    plt.savefig("../output/results/bifurcation_analysis/cenario" + scenario + "/itv=" + str(itv) + ".png")
    plt.show()
    plt.close()