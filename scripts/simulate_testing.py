from utils import Parameters, OutputData, write_testing_parameters, get_day_to_start_testing
from models import SeahirModel
import numpy as np
import os
import csv
import time

# Simulações para caso selvagem - salva tabelas de xI

def save_tables(tables, uf, num):
    filepath = os.path.join('..', 'output', 'testing', 'xI_table_'+uf+'_'+str(num)+'.csv')
    with open(filepath, 'w') as csv_file:
        spamwriter = csv.writer(csv_file)
        spamwriter.writerow(np.concatenate((['xI'], ['maxH'], ['mortos'], ['duracao'])))
        for i in range (0, 100):
            spamwriter.writerow(np.concatenate(([str(tables[i, 0])], [str(tables[i, 1])], [str(tables[i, 2])], [str(tables[i, 3])])))

# Definitions
I0 = 1
uf = 'SP'
I_to_start_testing = 10
parameters = Parameters(uf)
model = SeahirModel(uf)

I_set = [10, 100, 1000, 10000]
#Tabela para guardar, para cada I: maxH, mortos, duracao por xI
result_tables = np.zeros([4, 100, 4], dtype=np.float64)
age_strata = 16

for i, I_to_start_testing in enumerate(I_set):
    # First run - discover day with X infected
    parameters.reset()
    model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1, 1, 1], days=[1, 2, 3],
                   itv=[0, 0, 0], ifrs=parameters.ifrs, testing_parameters=False, verbose=True)
    # Discover day to start testing
    day_to_start_testing = get_day_to_start_testing(uf, I_to_start_testing, parameters)
    # print("Starting day "+str(day_to_start_testing))
    # Runs start testing on day X
    succesfull_xI = 0
    for j, xI in enumerate(np.arange(0, 1.0, 0.01)):
        print("Analyzing " + str(xI)+ " of "+str(I_to_start_testing))
        # write_testing_parameters(uf, xI=xI, xA=0)
        xI_array = np.full(age_strata, xI, dtype=np.float64)
        xA_array = np.zeros(age_strata, dtype=np.float64)
        model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1], days=[day_to_start_testing],
                   n=2, itv=[0], ifrs=parameters.ifrs, xI=xI_array, xA=xA_array, testing_parameters=False, verbose=False)
        outData = OutputData(uf, parameters)
        maxH = np.max(outData.H)
        deaths = outData.CD[-1]
        duration = np.max(np.argwhere(outData.I > 1))
        result_tables[i, j] = np.array([xI, maxH, deaths, duration])
    save_tables(result_tables[i], uf, I_to_start_testing)
        # if outData.I[-1] < 1:
        #     succesfull_xI = xI
        #     break

# print("xI minimo para fim da pandemia: " + str(succesfull_xI))
# plot_fit(uf, parameters, day_to_start_testing)
