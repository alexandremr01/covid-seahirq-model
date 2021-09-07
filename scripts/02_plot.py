import numpy as np
import matplotlib.pyplot as plt
import itertools

I0 = 1
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [0, 6, 8, 9]

for scenario, itv in list(itertools.product(scenarios, interventions)):
    y = np.loadtxt("../output/results/2_R0_xI_relation/cenario" + scenario + "/itv=" + str(itv) + ".csv", delimiter=',')
    R = y[:, 1]
    xI = y[:, 0]
    plt.plot(xI, R)
    plt.xlabel('xI')
    plt.ylabel('R0')
    name = ""

    if scenario == 'BR':
        name = "Brasil"
    elif scenario == 'GR':
        name = "Alemanha"
    elif scenario == 'NG':
        name = "Nigéria"
    elif scenario == 'US':
        name = "EUA"

    if itv == 0:
        plt.title(name+' - Sem intervenção')
    elif itv == 6:
        plt.title(name+' - Isolamento Vertical')
    elif itv == 8:
        plt.title(name+' - Distanciamento social')
    elif itv == 9:
        plt.title(name+' - Distanciamento social com restrições no trabalho')
        
    plt.savefig("../output/results/2_R0_xI_relation/cenario" + scenario + "/out_itv=" + str(itv) + ".png")
    # plt.show()
    plt.close()