from utils import OutputData, write_testing_parameters, get_day_to_start_testing, second_derivatives
import numpy as np
import os
import itertools
import pandas as pd
import matplotlib.pyplot as plt

I0 = 1
scenario = 'BR'
xi_set = [0.2, 0.5, 0.7, 0.9]
age_strata=16
default_figsize = [6.4, 4.8]

fig, axs = plt.subplots(4, 4, figsize=(4*default_figsize[0], 4*default_figsize[1]))

for i, xA in enumerate(xi_set):
    for j, xI in enumerate(xi_set):    # Run once without testing to get reference values
        x = np.array(list(range(400)))
        y = np.loadtxt("../output/results/03_death_day_relation/scenario"+scenario+"/xI=%.2f_xA=%.2f.csv" % (xI, xA), delimiter=',')
        y = y[:, 1]
        no_test_deaths = y[-1]

        y_xx = second_derivatives(y)
        print(y_xx)

        aux = np.argwhere(y_xx < -0.1)
        if len(aux)>1 :
            transition_day = np.min(aux)
        else:
            transition_day = 0

        y_np = np.array(y)
        y_np = y_np / no_test_deaths
        textstr = '\n'.join(
            [r'Óbitos sem testagem:%.0f' % (np.max(y)), 'Óbitos com testagem no D0: %.0f' % (np.min(y)), 'Óbitos evitáveis: %.0f' % (np.max(y) - np.min(y)),
            'Variação: {:.2%} '.format((np.max(y) - np.min(y)) / np.max(y)),
            'Dia de virada da 2a derivada: %d' % (transition_day)])
            #  'Infectados na virada=%.2e' % (outData.I[transition_day])])
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # axs[i, j].text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
        y_np[0] = np.min(y)
        # axs[i, j].axvline(x=transition_day, ymin=0, ymax=1.1 * np.max(y), color='black')

        if scenario == 'BR':
            name = "Brasil"
        elif scenario == 'GR':
            name = "Alemanha"
        elif scenario == 'NG':
            name = "Nigéria"
        elif scenario == 'US':
            name = "EUA"
            
        # y_np = (y_np-np.min(y_np))/(np.max(y_np)-np.min(y_np))
        # plt.scatter(x, y_np)
        # plt.scatter(x, I_norm)
        # plt.xlabel('Dia de inicio dos testes')
        # plt.title('Número de mortos ao final da epidemia com número de infectados normalizados\nPor dia de início das testagens com xI='+str(xI))
        # plt.savefig('../output/results/1_test_start/cenario'+scenario+'/Mortos por dia de inicio da testagem e infectados normalizados xI='+str(xI)+'.png')
        # # plt.show()
        # plt.close()
        # plt.axvline(x=transition_day, ymin=0, ymax= 1.1 * np.max(y), color='black')

        # axs[i, j].text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
        y[0] = np.min(y)
        axs[i, j].plot(x, y)
        red =  100*( no_test_deaths-y[0])/no_test_deaths
        # axs[i, j].xlabel('Dia de inicio dos testes')
        # axs[i, j].ylabel('Número de mortos ao final da epidemia')
        axs[i, j].set_title('xA=%s xI=%s\nRedução de %.2f%% de óbitos' % (str(xA), str(xI), red), {'fontsize': 20})


for ax in axs.flat:
    ax.set_xlabel('Dia de início dos testes', fontsize=18)
    ax.set_ylabel('Número de mortos ao final', fontsize=18)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

fig.tight_layout()
plt.savefig("../output/results/03_death_day_relation/result.png", dpi=300)

# axs[i, j].savefig('../output/results/1_test_start/cenario'+scenario+'/Mortos por dia de inicio da testagem xI='+str(xI)+'.png')
