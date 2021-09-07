import numpy as np
import matplotlib.pyplot as plt
import itertools
import matplotlib.cm as cm
import matplotlib.colors as mcolors

offset = mcolors.TwoSlopeNorm( vcenter=1.00, vmin=0.48, vmax=2.43)

I0 = 1
scenarios = ['BR', 'GR', 'NG', 'US']
age_strata = 16
interventions = [0, 6, 8, 9]
default_figsize = [6.4, 4.8]
fig, axs = plt.subplots(4, 4, figsize=(4*default_figsize[0], 4*default_figsize[1]))
cmap=cm.get_cmap('bwr')

Rmins=[]
epsilon = 1e-3

for i, scenario in enumerate(scenarios):
    for j, itv in enumerate(interventions):
        efficient_points_x = []
        efficient_points_y=[]

        y = np.loadtxt("../output/results/06_R0_xi_xa/cenario" + scenario + "/itv=" + str(itv) + ".csv", delimiter=',')
        R_raw = y[:, 2]
        xI = y[:, 0]
        xA = y[:, 1]
        m = 100

        R = np.zeros((m, m), dtype=np.float64)
        for p in range(m):
            for q in range(m):
                R[p, q] = R_raw[p*100 + q]
                
                if abs(R_raw[p*100 + q] - 1) < epsilon:
                    efficient_points_x.append(xA[p*100 + q])
                    efficient_points_y.append(xI[p*100 + q])

        Rmins.append(np.min(R))

        pcm = axs[i, j].imshow(R, extent=[np.min(xA), np.max(xA), np.min(xI), np.max(xI)], origin='lower', interpolation='none', cmap=cmap, norm=offset)
        axs[i, j].plot(efficient_points_x, efficient_points_y, color='black')
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
            itvname = ' - Restrições no trabalho'
        axs[i, j].set_title(name+itvname, {'fontsize': 20})

for ax in axs.flat:
    ax.set_xlabel('xA', fontsize=18)
    ax.set_ylabel('xI', fontsize=18)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

fig.tight_layout()

print(np.min(Rmins))
im=cm.ScalarMappable(cmap=cmap, norm=offset)
cb=  fig.colorbar(im, ax=axs.ravel().tolist())

cb.set_label("R0", fontsize=15)
plt.savefig("../output/results/06_R0_xi_xa/result.png", dpi=300)