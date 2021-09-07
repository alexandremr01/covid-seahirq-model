# coding: utf-8
import numpy as np
import csv
import os
import subprocess
from math import log

compartments = 12
num_intervention_phases = 3
age_strata = 16
t_days = 400

def second_derivative(x, t):
    if t>0 and t<len(x)-1:
        return x[t-1] - 2*x[t] + x[t+1]
    else:
        return 0

def second_derivatives(x):
    b = np.zeros(x.shape, dtype=np.float64)
    for t in range(len(x)):
        b[t] = second_derivative(x, t)
    return b

def first_derivative(x, t):
    if t>0 and t<len(x)-1:
        return 0.5*(x[t+1] - x[t-1])
    elif t==0:
        return x[t+1]-x[t]
    elif t==len(x)-1:
        return x[t]-x[t-1]

def first_derivatives(x):
    b = np.zeros(x.shape, dtype=np.float64)
    for t in range(len(x)):
        b[t] = first_derivative(x, t)
    return b

def doubling_time(y, t, window = 10):
    y1 = y[t]
    y2 = y[t+window]
    return window * log(2) / log(y2 / y1)

def get_day_to_start_testing(uf, i0, parameters):
    outData = OutputData(uf, parameters)
    return np.argmax(outData.I > i0) + 2

def write_testing_parameters(uf, xI, xA):
    testing_parametres_file = os.path.join('..', '..', 'input', 'cenarios', 'cenario'+uf, 'testing_parameters.csv')
    with open(testing_parametres_file, 'w') as csv_file:
        spamwriter = csv.writer(csv_file)

        spamwriter.writerow(['VARIAVEL', 'FAIXA_1', 'FAIXA_2', 'FAIXA_3', 'FAIXA_4', 'FAIXA_5', 'FAIXA_6', 'FAIXA_7', 'FAIXA_8',
                'FAIXA_9', 'FAIXA_10', 'FAIXA_11', 'FAIXA_12', 'FAIXA_13', 'FAIXA_14', 'FAIXA_15', 'FAIXA_16'])
        spamwriter.writerow(np.concatenate((['xI0'], 16*[str('0')])))
        spamwriter.writerow(np.concatenate((['xA0'], 16*[str('0')])))
        spamwriter.writerow(np.concatenate((['xI1'], 16*[str(xI)])))
        spamwriter.writerow(np.concatenate((['xA1'], 16*[str(xA)])))

class OptimizedParameters:
    """
    Lê os parâmetros salvos.
    """

    def __init__(self, uf, num_intervention_phases=3, tuned=True):
        self.uf = uf
        if tuned:
            cenario = '../../input/cenarios/cenario' + uf + '/optimized_parameters_tuned.csv'
            parameters = []
            with open(cenario, "r") as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                for row in spamreader:
                    if len(row) > 0:
                        parameters.append((row[1]))
            self.r_0 = float(parameters[0])
            self.tau = float(parameters[1])
            self.zeta = float(parameters[2])
            self.D0 = float(parameters[3])
            n = num_intervention_phases
            self.g = [float(x) for x in parameters[4: n + 4]]
            durations = [float(x) for x in parameters[n + 4 : 2*n+4]]
            self.days = np.add.accumulate(durations)
            self.ifrs = [float(x) for x in parameters[2 * n + 4:]]
        self.gamma_H = self.read_gamma_h()
        # print(parameters)

    def read_gamma_h(self):
        param_file = '../../input/cenarios/cenario' + self.uf + '/parameters.csv'
        with open(param_file, "r") as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                if len(row) > 0 and row[0] == "GAMA_H":
                    gamma_h = np.array(row[1:], dtype=np.float)
        # print(gamma_h)
        return gamma_h.reshape(1, -1)

    def reset(self, g=[1, 1, 1], days=[1, 2, 3], itv=[0, 0, 0]):
        self.g = g
        self.days = days
        self.itv = itv


class OutputData:
    """
    Lê o output do modelo.
    """

    def __init__(self, uf, param=None):
        self.s = np.zeros([t_days, age_strata], dtype=np.float64)
        self.e = np.zeros([t_days, age_strata], dtype=np.float64)
        self.y = np.zeros([t_days, age_strata], dtype=np.float64)
        self.r = np.zeros([t_days, age_strata], dtype=np.float64)
        self.n = np.zeros([t_days, age_strata], dtype=np.float64)
        self.a = np.zeros([t_days, age_strata], dtype=np.float64)
        self.c = np.zeros([t_days, age_strata], dtype=np.float64)
        self.d = np.zeros([t_days, age_strata], dtype=np.float64)
        self.h = np.zeros([t_days, age_strata], dtype=np.float64)
        self.l = np.zeros([t_days, age_strata], dtype=np.float64)
        self.ri = np.zeros([t_days, age_strata], dtype=np.float64)
        self.qi = np.zeros([t_days, age_strata], dtype=np.float64)
        self.qa = np.zeros([t_days, age_strata], dtype=np.float64)
        self.t_array = np.array(t_days, dtype=np.float64)

        self.t_array = np.linspace(0, t_days - 1, t_days)

        output_file = '../../output/result_data_' + uf + '.csv'

        with open(output_file, "r") as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            next(spamreader, None)
            j = 0
            for row in spamreader:
                for i in range(age_strata):
                    self.s[j, i] = row[compartments * (i + 1)]
                    self.e[j, i] = row[compartments * (i + 1) + 1]
                    self.y[j, i] = row[compartments * (i + 1) + 2]
                    self.r[j, i] = row[compartments * (i + 1) + 3]
                    self.n[j, i] = row[compartments * (i + 1) + 4]
                    self.a[j, i] = row[compartments * (i + 1) + 5]
                    self.c[j, i] = row[compartments * (i + 1) + 6]
                    self.d[j, i] = row[compartments * (i + 1) + 7]
                    self.h[j, i] = row[compartments * (i + 1) + 8]
                    self.l[j, i] = row[compartments * (i + 1) + 9]
                    self.qi[j, i] = row[compartments * (i + 1) + 10]
                    self.qa[j, i] = row[compartments * (i + 1) + 11]
                for ii in range(age_strata):
                    self.ri[j, ii] = row[compartments * (age_strata + 1) + ii + 1]
                j = j + 1
        self.cd = self.c + self.d
        self.S = np.sum(self.s, axis=1)
        self.E = np.sum(self.e, axis=1)
        self.I = np.sum(self.y, axis=1)
        self.R = np.sum(self.r, axis=1)
        self.N = np.sum(self.n, axis=1)
        self.A = np.sum(self.a, axis=1)
        self.C = np.sum(self.c, axis=1)
        self.D = np.sum(self.d, axis=1)
        self.H = np.sum(self.h, axis=1)
        self.Qi = np.sum(self.qi, axis=1)
        self.Qa = np.sum(self.qa, axis=1)
        # print(param.gamma_h.shape)
        # print(self.y.shape)
        if param != None:
            self.hac = np.add.accumulate(np.multiply(param.gamma_H, self.y), axis=0)
            self.HAc = np.sum(self.hac, axis=1)
        self.L = np.sum(self.l, axis=1)
        self.RI = np.sum(self.ri, axis=1)
        self.CD = self.C + self.D
        # print(self.CD)
        self.Ac = self.C + self.RI + self.I + self.D
