from utils import Parameters, OutputData, write_testing_parameters, get_day_to_start_testing
from models import SeahirModel
import numpy as np
from plot_test_simulations import plot_fit
import os
import csv

# Simulações para caso selvagem - salva tabelas de xI

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

I_set = [10, 100, 1000, 10000, 100000, 1000000, 3000000]
xI = 0.5

result_tables = np.zeros([len(I_set), 4], dtype=np.float64)

for i, I_to_start_testing in enumerate(I_set):
    # First run - discover day with X infected
    parameters.reset()
    model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1, 1, 1], days=[1, 2, 3],
                   itv=[0, 0, 0], ifrs=parameters.ifrs, testing_parameters=False, verbose=False)
    outData = OutputData(uf, parameters)
    maxI = np.argmax(outData.I)
    print("I max: ", maxI)
    # Discover day to start testing
    day_to_start_testing = get_day_to_start_testing(uf, I_to_start_testing, parameters)
    # Runs start testing on day X
    write_testing_parameters(uf, xI=xI, xA=0)
    model.run_ifrs(I0=I0, R0=2.0, zeta=parameters.zeta, tau=parameters.tau, g=[1], days=[day_to_start_testing],
               n=2, itv=[0], ifrs=parameters.ifrs, testing_parameters=True, verbose=False)
    outData = OutputData(uf, parameters)
    maxH = np.max(outData.H)
    deaths = outData.CD[-1]
    duration = np.max(np.argwhere(outData.I > 1))
    result_tables[i] = np.array([I_to_start_testing, maxH, deaths, duration])
save_tables(result_tables, uf, len(I_set))
