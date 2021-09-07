import numpy as np
import matplotlib.pyplot as plt
import itertools

I0 = 1
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [0, 6, 8, 9]
default_figsize = [6.4, 4.8]
fig, axs = plt.subplots(4, 4, figsize=(4*default_figsize[0], 4*default_figsize[1]))

for i, scenario in enumerate(scenarios):
    for j, itv in enumerate(interventions):
        y = np.loadtxt("../output/results/2_R0_xI_relation/cenario" + scenario + "/itv=" + str(itv) + ".csv", delimiter=',')
        R = y[:, 1]
        xI = y[:, 0]
        axs[i, j].plot(xI, R)

        if scenario == 'BR':
            name = "Brasil"
        elif scenario == 'GR':
            name = "Alemanha"
        elif scenario == 'NG':
            name = "Nigéria"
        elif scenario == 'US':
            name = "EUA"

        if itv == 0:
            itvname = ' - Sem intervenção'
        elif itv == 6:
            itvname = ' - Isolamento Vertical'
        elif itv == 8:
            itvname = ' - Distanciamento social'
        elif itv == 9:
            itvname = ' - Com restrições no trabalho'
        axs[i, j].set_title(name+itvname, {'fontsize': 20})

for ax in axs.flat:
    ax.set_xlabel('xI', fontsize=18)
    ax.set_ylabel('R0', fontsize=18)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

fig.tight_layout()
plt.savefig("../output/results/2_R0_xI_relation/result.png", dpi=300)