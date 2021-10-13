import numpy as np
import matplotlib.pyplot as plt
import itertools
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy import interpolate

I0 = 1
scenarios = ['BR', 'GR', 'US', 'NG']
interventions = [0, 6, 8, 9]
# scenarios = [ 'US']
# interventions = [0]
age_strata = 16
default_figsize = [6.4, 4.8]
fig, axs = plt.subplots(4, 4, figsize=(4*default_figsize[0], 4*default_figsize[1]))
cmap=cm.get_cmap('Reds')

Rmins=[]
epsilon = 1e-3
xa_values = np.linspace(0, 0.98, 50)
step = 0.98/49

from scipy import interpolate
def interpolate_missing_pixels(
        image: np.ndarray,
        mask: np.ndarray,
        method: str = 'nearest',
        fill_value: int = 0
):
    """
    :param image: a 2D image
    :param mask: a 2D boolean image, True indicates missing values
    :param method: interpolation method, one of
        'nearest', 'linear', 'cubic'.
    :param fill_value: which value to use for filling up data outside the
        convex hull of known pixel values.
        Default is 0, Has no effect for 'nearest'.
    :return: the image with missing values interpolated
    """
    h, w = image.shape[:2]
    xx, yy = np.meshgrid(np.arange(w), np.arange(h))
    known_x = xx[~mask]
    known_y = yy[~mask]
    known_v = image[~mask]
    missing_x = xx[mask]
    missing_y = yy[mask]
    interp_values = interpolate.griddata(
        (known_x, known_y), known_v, (missing_x, missing_y),
        method=method, fill_value=fill_value
    )
    interp_image = image.copy()
    interp_image[missing_y, missing_x] = interp_values
    return interp_image

reverse_map = {}
for i in range(len(xa_values)):
    reverse_map[xa_values[i]] = i

for i, scenario in enumerate(scenarios):
    for j, itv in enumerate(interventions):
        efficient_points_x = []
        efficient_points_y = []
        print(itv)
        y = np.loadtxt("../output/results/02_critical_attack/scenario" + scenario + "/itv=" + str(itv) + ".csv", delimiter=',')
        critical_attack_raw = y[:, 2]
        xI = y[:, 0]
        xA = y[:, 1]
        m = 50
        critical_attack = np.zeros((m, m), dtype=np.float64)
        for k in range(len(critical_attack_raw)):
            critical_attack[round(xI[k]/step), round(xA[k]/step)] = critical_attack_raw[k]
        critical_attack[critical_attack == 0.0] = np.nan
        # critical_attack = interpolate_missing_pixels(critical_attack, np.isnan(critical_attack), method='nearest')

        pcm = axs[i, j].imshow(critical_attack, extent=[np.min(xA), np.max(xA), np.min(xI), np.max(xI)], origin='lower', interpolation='none', cmap=cmap)
        cm = fig.colorbar(pcm, ax=axs[i, j], location='right', anchor=(0, 0.3))
      


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

# im=cm.ScalarMappable(cmap=cmap)
# cb=  fig.colorbar(im, ax=axs.ravel().tolist())

# cb.set_label("Susceptibilidade Crítica", fontsize=15)
plt.savefig("../output/results/02_critical_attack/result.png", dpi=300)