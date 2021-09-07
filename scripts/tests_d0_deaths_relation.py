from utils import Parameters, OutputData, write_testing_parameters, get_day_to_start_testing
from models import SeahirModel
import numpy as np
from plot_test_simulations import plot_fit
import os
import csv
import matplotlib.pyplot as plt
# Correlaciona dia de inicio dos testes, com xI fixo, e número de mortos

def save_tables(tables, uf, qnt):
    filepath = os.path.join('..', 'output', 'testing', 'xI_I0_table_'+uf+'.csv')
    with open(filepath, 'w') as csv_file:
        spamwriter = csv.writer(csv_file)
        spamwriter.writerow(np.concatenate((['I1'], ['maxH'], ['mortos'], ['duracao'])))
        for i in range (0, qnt):
            spamwriter.writerow(np.concatenate(([str(tables[i, 0])], [str(tables[i, 1])], [str(tables[i, 2])], [str(tables[i, 3])])))

# Definitions
I0 = 1
uf = 'SP'
I_to_start_testing = 10
parameters = Parameters(uf)
model = SeahirModel(uf)
plot_infections = True
I_set = [10, 100, 1000, 10000, 100000, 1000000, 3000000]
xi_set = [0.2, 0.5, 0.7, 0.9]
for xI in xi_set:
    write_testing_parameters(uf, xI=xI, xA=0)

    model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1], days=[0],
                   n=2, itv=[0], ifrs=parameters.ifrs, testing_parameters=False, verbose=False)
    outData = OutputData(uf, parameters)
    I = outData.I[0:400]/np.max(outData.I[0:400])
    x = list(range(400))
    y = []

    age_strata = 16
    xI_array = np.full(age_strata, xI, dtype=np.float64)
    xA_array = np.zeros(age_strata, dtype=np.float64)
    for d in range(400):
        print("Running day ", d)
        # Runs start testing on day d
        # Com testing parameters
        # model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1], days=[d],
        #            n=2, itv=[0], ifrs=parameters.ifrs, testing_parameters=True, verbose=False)
        model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1], days=[d],
                   n=2, itv=[0], ifrs=parameters.ifrs, testing_parameters=False, xI=xI_array, xA=xA_array, verbose=False)
        outData = OutputData(uf, parameters)
        maxH = np.max(outData.H)
        deaths = outData.CD[-1]
        duration = np.max(np.argwhere(outData.I > 1))
        y.append(deaths)

    np.savetxt("../output/testing/per_day/data_xI="+str(xI)+".csv", y, delimiter=",")

    # y_np = np.array(y)
    # y_np = (y_np-np.min(y_np))/(np.max(y_np)-np.min(y_np))
    # plt.scatter(x, y_np)
    # plt.scatter(x, I)
    # plt.xlabel('Dia de inicio dos testes')
    # plt.title('Número de mortos ao final da epidemia com número de infectados normalizados\nPor dia de início das testagens com xI='+str(xI))
    # plt.savefig('../output/testing/per_day/Mortos por dia de inicio da testagem e infectados normalizados xI='+str(xI)+'.png')
    # #plt.show()
    # plt.close()
    #
    # plt.scatter(x, y)
    # plt.xlabel('Dia de inicio dos testes')
    # plt.ylabel('Número de mortos ao final da epidemia')
    # plt.title('Número de mortos ao final da epidemia\nxI='+str(xI))
    # plt.savefig('../output/testing/per_day/Mortos por dia de inicio da testagem xI='+str(xI)+'.png')
    #plt.show()
    # plt.show()