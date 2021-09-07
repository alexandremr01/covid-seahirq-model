from utils import OutputData, write_testing_parameters, get_day_to_start_testing, second_derivatives
from models import SeahirModel, EpidemicPhases
import numpy as np
from model.nextgen import calculateR0, R0FromBetaGama, R0FromxIxA
import os
import itertools
from model.scen_gen import scen_gen, TESTING_PARAMETERS_DIRECT
import pandas as pd
from model.parameters import Parameters
import matplotlib.pyplot as plt

I0 = 1
scenarios = ['GR']
xi_set = [0.2, 0.5, 0.7, 0.9]
age_strata=16

for uf in scenarios:
    # Run once without testing to get reference values
    model = SeahirModel(uf)
    phases = EpidemicPhases(g0=1, itv0=0)
    # phases.add_phase(g=1, start_day=0, itv_level=0)
    xI_array = np.full(age_strata, 0, dtype=np.float64)
    xA_array = np.zeros(age_strata, dtype=np.float64)

    model.run(I0=I0, R0=None, phases=phases, verbose=False, attack=1.0, xI=xI_array, xA=xA_array, testing_parameters = TESTING_PARAMETERS_DIRECT)
    outData = OutputData(uf)
    I_norm = outData.I[0:400] / np.max(outData.I[0:400])

    for xI in xi_set:
        x = np.array(list(range(400)))
        y = np.loadtxt("../output/results/1_test_start/cenario"+uf+"/data_xI="+str(xI)+".csv", delimiter=',')
        y = y[:, 1]

        y_xx = second_derivatives(y)
        transition_day = np.min(np.argwhere(y_xx < -0.01))

        y_np = np.array(y)
        textstr = '\n'.join(
            [r'Óbitos sem testagem:%.0f' % (np.max(y)), 'Óbitos com testagem no D0=%.0f' % (np.min(y)), 'Óbitos evitáveis=%.0f' % (np.max(y) - np.min(y)),
             'Variação={:.2%} '.format((np.max(y) - np.min(y)) / np.min(y)),
             'Dia de virada da 2a derivada=%d' % (transition_day)])
            #  'Infectados na virada=%.2e' % (outData.I[transition_day])])
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
        y_np[0] = np.min(y)
        plt.axvline(x=transition_day, ymin=0, ymax=1.1 * np.max(y), color='black')

        if uf == 'BR':
            name = "Brasil"
        elif uf == 'GR':
            name = "Alemanha"
        elif uf == 'NG':
            name = "Nigéria"
        elif uf == 'US':
            name = "EUA"
            
        y_np = (y_np-np.min(y_np))/(np.max(y_np)-np.min(y_np))
        plt.scatter(x, y_np)
        plt.scatter(x, I_norm)
        plt.xlabel('Dia de inicio dos testes')
        plt.title('Número de mortos ao final da epidemia com número de infectados normalizados\nPor dia de início das testagens com xI='+str(xI))
        plt.savefig('../output/results/1_test_start/cenario'+uf+'/Mortos por dia de inicio da testagem e infectados normalizados xI='+str(xI)+'.png')
        # plt.show()
        plt.close()
        plt.axvline(x=transition_day, ymin=0, ymax= 1.1 * np.max(y), color='black')

        plt.text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
        y[0] = np.min(y)
        plt.scatter(x, y)
        plt.xlabel('Dia de inicio dos testes')
        plt.ylabel('Número de mortos ao final da epidemia')
        plt.title('Número de mortos ao final da epidemia\n%s xI=%s' % (name, str(xI)))
        plt.savefig('../output/results/1_test_start/cenario'+uf+'/Mortos por dia de inicio da testagem xI='+str(xI)+'.png')
        # plt.show()
        plt.close()