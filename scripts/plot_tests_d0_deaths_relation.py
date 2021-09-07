import numpy as np
import matplotlib.pyplot as plt
from utils import second_derivatives, Parameters, OutputData
xi_set = [0.2, 0.5, 0.7, 0.9]
from models import SeahirModel

uf = 'SP'
I_to_start_testing = 10
parameters = Parameters(uf)
model = SeahirModel(uf)

model.run_ifrs(I0=1, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1], days=[0],
                   n=2, itv=[0], ifrs=parameters.ifrs, testing_parameters=False, verbose=False)
outData = OutputData(uf, parameters)
I = outData.I[0:400] / np.max(outData.I[0:400])
for xI in xi_set:
    x = np.array(list(range(400)))
    y = np.loadtxt("../output/testing/per_day/data_xI="+str(xI)+".csv")

    y_xx = second_derivatives(y)
    transition_day = np.min(np.argwhere(y_xx < -0.01))


    y_np = np.array(y)
    textstr = '\n'.join(
        [r'$Max=%.2e$' % (np.max(y)), '$Min=%.2e$' % (np.min(y)), '$Delta=%.2e$' % (np.max(y) - np.min(y)),
         'Variação={:.2%} '.format((np.max(y) - np.min(y)) / np.min(y)),
         'Dia de virada da 2a derivada=%d' % (transition_day),
         'Infectados na virada=%.2e' % (outData.I[transition_day])])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
    y_np[0] = np.min(y)
    plt.axvline(x=transition_day, ymin=0, ymax=1.1 * np.max(y), color='black')
    y_np = (y_np-np.min(y_np))/(np.max(y_np)-np.min(y_np))
    plt.scatter(x, y_np)
    plt.scatter(x, I)
    plt.xlabel('Dia de inicio dos testes')
    plt.title('Número de mortos ao final da epidemia com número de infectados normalizados\nPor dia de início das testagens com xI='+str(xI))
    plt.savefig('../output/testing/per_day/Mortos por dia de inicio da testagem e infectados normalizados xI='+str(xI)+'.png')
    # plt.show()
    plt.close()
    plt.axvline(x=transition_day, ymin=0, ymax= 1.1 * np.max(y), color='black')

    plt.text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
    y[0] = np.min(y)
    plt.scatter(x, y)
    plt.xlabel('Dia de inicio dos testes')
    plt.ylabel('Número de mortos ao final da epidemia')
    plt.title('Número de mortos ao final da epidemia\nxI='+str(xI))
    plt.savefig('../output/testing/per_day/Mortos por dia de inicio da testagem xI='+str(xI)+'.png')
    # plt.show()
    plt.close()