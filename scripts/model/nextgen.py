# -*- coding: utf-8 -*-
import sympy as sym
import numpy as np
import sys
from scipy import linalg
import pandas as pd

# Infected columns
# Exposed: de 0*12 a 1*12
# X: de x*12 a (x+1)*12
E = 0
A = 1
I = 2
Qi = 3
Qa = 4
H = 5

age_strata = 16

def R0FromBetaGama(scenario, changeDir = True, verbose = False):
    if changeDir:
        input_folder = '../input/cenarios/cenario' + scenario + '/'
    else:
        input_folder = ''
    df = pd.read_csv(input_folder + 'beta_gama.csv', sep=',')
    xI = df.values[0, 17:17 + age_strata]
    xA = df.values[0, 17 + age_strata:17 + 2 * age_strata]
    betas = df.values[0:age_strata, 17 + 2 * age_strata:17 + 3 * age_strata]
    return calculateR0(scenario, betas, xI, xA, changeDir, verbose)


def R0FromxIxA(scenario, xI, xA, attack=1.0, changeDir = True, verbose = False, nextGenScenario=None):
    return calculateR0(scenario, None, xI, xA, changeDir, verbose, nextGenScenario=nextGenScenario, attack=attack)

class NextGenScenario:
    def __init__(self, scenario, changeDir):
        if changeDir:
            input_folder = '../input/cenarios/cenario' + scenario + '/'
        else:
            input_folder = ''
        # Read parameters
        df = pd.read_csv(input_folder + 'parameters.csv', sep=',')

        self.mu_eq = df.values[1, 1:]
        self.mu_cov = df.values[2, 1:]
        self.theta = df.values[3, 1:]
        self.alphas = df.values[4, 1:]
        self.ksi = df.values[5, 1:]
        self.rho = df.values[6, 1:]
        self.phi = df.values[7, 1:]
        self.etas = df.values[8, 1:]
        self.a = df.values[9, 1:]
        self.gamma_H = df.values[10, 1:]
        self.gamma_HR = df.values[11, 1:]
        self.gamma_RI = df.values[12, 1:]
        self.gamma_RA = df.values[13, 1:]
        self.gamma_RQI = df.values[14, 1:]
        self.gamma_RQA = df.values[15, 1:]
        self.gamma_HQI = df.values[16, 1:]

        #Read demography
        df = pd.read_csv(input_folder + 'demographic_data.csv', sep=',')
        pop_strata = df.values[0, 1:]
        self.pop_strata = np.divide(pop_strata,np.sum(pop_strata))

        df = pd.read_csv(input_folder + 'beta_gama.csv', sep=',')
        self.betas = df.values[17:17+age_strata, 17 + 2 * age_strata:17 + 3 * age_strata]



def calculateR0(scenario, contact_matrix, xI, xA, changeDir = True, verbose = False, nextGenScenario=None, attack=1.0):
    if nextGenScenario is None:
        nextGenScenario = NextGenScenario(scenario, changeDir)
    else: 
        contact_matrix = attack*nextGenScenario.betas

    mu_eq = nextGenScenario.mu_eq
    mu_cov = nextGenScenario.mu_cov
    theta = nextGenScenario.theta
    alphas = nextGenScenario.alphas
    ksi = nextGenScenario.ksi
    rho = nextGenScenario.rho
    phi = nextGenScenario.phi
    etas = nextGenScenario.etas
    a = nextGenScenario.a
    gamma_H = nextGenScenario.gamma_H
    gamma_HR = nextGenScenario.gamma_HR
    gamma_RI = nextGenScenario.gamma_RI
    gamma_RA = nextGenScenario.gamma_RA
    gamma_RQI = nextGenScenario.gamma_RQI
    gamma_RQA = nextGenScenario.gamma_RQA
    gamma_HQI = nextGenScenario.gamma_HQI
    pop_strata = nextGenScenario.pop_strata

    # Read beta gamma
    betas = contact_matrix

    gamma_QI = np.multiply(xI, gamma_H+gamma_RI) / (1-xI)
    gamma_QA = np.multiply(xA, gamma_RA) / (1-xA)

    # A unica entrada de novas infecções é em E, vindo de A, I ou E
    F = np.zeros((6*age_strata, 6*age_strata))

    # Some matrix manipulation
    B_1 = np.multiply(ksi, np.multiply(betas, np.tile(pop_strata.reshape(age_strata, 1), (1, age_strata))))
    B_2 = np.multiply(alphas, np.multiply(betas, np.tile(pop_strata.reshape(age_strata, 1), (1, age_strata))))
    B_3 = np.multiply(betas, np.tile(pop_strata.reshape(age_strata, 1), (1, age_strata)))
    F[E * age_strata: (E + 1) * age_strata, E * age_strata: (E + 1) * age_strata] = B_1
    F[E * age_strata: (E + 1) * age_strata, A * age_strata: (A + 1) * age_strata] = B_2
    F[E * age_strata: (E + 1) * age_strata, I * age_strata: (I + 1) * age_strata] = B_3

    # Transições entre compartimentos de infectados
    V = np.zeros((6*age_strata, 6*age_strata))
    # Exposto -> Exposto
    for i in range(age_strata):
        j = E + i
        V[j, j] = mu_eq[i] + a[i]
    # Assintomático
    for i in range(age_strata):
        j = A*age_strata + i
        V[j, j] = mu_eq[i] + gamma_RA[i] + gamma_QA[i]
        V[j, E*age_strata+i] = -a[i] * (1-rho[i])
    # Infectado sintomático
    for i in range(age_strata):
        j = I*age_strata + i
        V[j, j] = mu_eq[i] + theta[i]*mu_cov[i] + gamma_H[i] + gamma_RI[i] + gamma_QI[i]
        V[j, E*age_strata+i] = -a[i]*rho[i]
    # Quarentenado sintomático
    for i in range(age_strata):
        j = Qi*age_strata + i
        V[j, j] = mu_eq[i] + gamma_HQI[i] + gamma_RQI[i]
        V[j, I*age_strata+i] = -gamma_QI[i]
    # Quarentenado assintomático
    for i in range(age_strata):
        j = Qa*age_strata + i
        V[j, j] = mu_eq[i] + gamma_RQA[i]
        V[j, A*age_strata+i] = -gamma_QA[i]
    # Hospitalizado
    for i in range(age_strata):
        j = H*age_strata + i
        V[j, j] = mu_eq[i] + gamma_HR[i] + mu_cov[i]
        V[j, I*age_strata+i] = -gamma_H[i]
        V[j, Qi*age_strata+i] = -gamma_HQI[i]

    # print(V)
    v_1 = np.linalg.inv(V)
    np.set_printoptions(threshold=sys.maxsize)
    # np.savetxt("V.csv", V, delimiter=",")
    # np.savetxt("F.csv", F, delimiter=",")
    next_gen_matrix = np.matmul(F, v_1)
    w, _ = linalg.eig(next_gen_matrix)
    R0 = np.max(np.abs(w))
    if verbose:
        print("R0 obtido: "+str(R0))
    return R0