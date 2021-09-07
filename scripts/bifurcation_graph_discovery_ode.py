from functools import partial 
from collections import defaultdict 
import numpy as np # Numerical computing library
import matplotlib.pyplot as plt # Plotting library
import scipy.integrate #Integration library
from mpl_toolkits.mplot3d import axes3d #Used for the 3d bifurcation plot
import matplotlib.patches as mpatches #used to write custom legends
from models import SeahirModel, EpidemicPhases
from model.scen_gen import scen_gen, TESTING_PARAMETERS_FROM_EPIDEMIOLOGY
import pandas as pd

scenario = 'BR'
age_strata = 16
model = SeahirModel(scenario)
input_folder = '../input/cenarios/cenarioBR/'


class Func:

    def __init__(self, R0):
        """ ODE system modeling Gardner's bistable cellular switch
        Args:
            y (array): (concentration of u, concentration of v)
            t (float): Time
            alpha (float): maximum rate of repressor synthesis 
            beta (float): degree of cooperative behavior.
        Return: dy/dt
        """
        print(R0)
        # 1. dado um R0, gera todos os parâmetros
        model.run(I0=1, R0=R0, phases=EpidemicPhases(g0=1, itv0=0), verbose=False, attack=None, testing_parameters = TESTING_PARAMETERS_FROM_EPIDEMIOLOGY)
        # 2. com os parâmetros, monta as EDOs
        # Read parameters
        df = pd.read_csv(input_folder + 'beta_gama.csv', sep=',')
        self.xI = df.values[0, 17:17 + age_strata]
        self.xA = df.values[0, 17 + age_strata:17 + 2 * age_strata]
        self.betas = df.values[0:age_strata, 17 + 2 * age_strata:17 + 3 * age_strata]

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
        self.pop_strata = df.values[0, 1:]
        self.natalidade = df.values[2, 1:]
        self.pop_strata = np.divide(self.pop_strata,np.sum(self.pop_strata))

        self.gamma_QI = np.multiply( self.xI, np.divide((self.gamma_H + self.gamma_RI), (1 - self.xI)) )
        self.gamma_QA = np.multiply( self.xA, np.divide(self.gamma_RA, (1 - self.xA) ))

    def get_f(self):
        def f(y, t=0):
            S = []
            E = []
            I = []
            R = []
            N = []
            A = []
            C = []
            D = []
            H = []
            L = []
            Qi = []
            Qa = []
            for i in range(age_strata):
                S.append(y[i + 0])
                E.append(y[i + 1])
                I.append(y[i + 2])
                R.append(y[i + 3])
                N.append(y[i + 4])
                A.append(y[i + 5])
                C.append(y[i + 6])
                D.append(y[i + 7])
                H.append(y[i + 8])
                L.append(y[i + 9])
                Qi.append(y[i + 10])
                Qa.append(y[i + 11])


            mu_eq = self.mu_eq
            mu_cov = self.mu_cov
            theta = self.theta
            alphas = self.alphas
            ksi = self.ksi
            rho = self.rho
            phi = self.phi
            etas = self.etas
            a = self.a
            gamma_H = self.gamma_H
            gamma_HR = self.gamma_HR
            gamma_RI = self.gamma_RI
            gamma_RA = self.gamma_RA
            gamma_RQI = self.gamma_RQI
            gamma_RQA = self.gamma_RQA
            gamma_HQI = self.gamma_HQI
            
            betas = self.betas
            pop_strata = self.pop_strata
            natalidade = self.natalidade
            pop_strata = self.pop_strata

            gamma_QI = self.gamma_QI
            gamma_QA = self.gamma_QA



            dS = []
            dE = []
            dA = []
            dI = []
            dQi = []
            dQa = []
            dH = []
            dR = []
            dRi = []
            dC = []
            dL = []
            dD = []
            dN = []
            for k in range(age_strata):
                population = np.sum(N)
                beta_row_sum = 0
                for j in range(age_strata):
                    beta_row_sum += betas[i, j] * S[k] * I[j] / population
                alpha_beta_row_sum = 0
                for j in range(age_strata):
                    alpha_beta_row_sum += alphas[j] * betas[i, j] * S[k] * I[j] / population
                exp_beta_row_sum = 0
                for j in range(age_strata):
                    exp_beta_row_sum += ksi[j] * betas[i, j] * S[k] * I[j] / population

                dS.append(natalidade[k]*N[k] - mu_eq[k] * S[k] - alpha_beta_row_sum - beta_row_sum - exp_beta_row_sum)
                dE.append( beta_row_sum + alpha_beta_row_sum + exp_beta_row_sum  - (mu_eq[k] + a[k]) * E[k])
                dA.append( a[k] * (1.0 - rho[k]) * E[k] - (mu_eq[k] + gamma_RA[k] + gamma_QA[k]) * A[k])
                dI.append( a[k] * rho[k] * E[k] - (theta[k] * mu_cov[k] + mu_eq[k] + gamma_H[k] + gamma_RI[k] + gamma_QI[k]) * I[k])
                dQi.append( gamma_QI[k] * I[k] - (gamma_HQI[k] + mu_eq[k] + gamma_RQI[k]) * Qi[k])
                dQa.append( gamma_QA[k] * A[k] - (gamma_RQA[k] + mu_eq[k]) * Qa[k])
                dH.append( gamma_H[k] * I[k] + gamma_HQI[k] * Qi[k] - (gamma_HR[k] + mu_eq[k] + mu_cov[k]) * H[k])
                dR.append( gamma_RI[k] * I[k] + gamma_RA[k] * A[k] + gamma_RQA[k] * Qa[k] + gamma_RQI[k] * Qi[k] + gamma_HR[k] * H[k] - mu_eq[k] * R[k])
                dC.append( mu_cov[k] * H[k])
                dD.append( theta[k] * mu_cov[k] *  I[k])
                dN.append( natalidade[k] * N[k] - mu_eq[k] * N[k] - mu_cov[k] * H[k] - theta[k] * mu_cov[k] * I[k])
                dL.append(0)

            return np.concatenate([dS, dE, dA, dI, dQi, dQa, dH, dR, dRi, dC, dD, dN, dL])
        return f


def findroot(func, init): 
    """ Find root of equation function(x)=0
    Args:
        - the system (function),
        - the initial values (type list or np.array)

    return: correct equilibrium (type np.array) 
            if the numerical method converge or return nan
    """
    sol, info, convergence, sms = scipy.optimize.fsolve(func, init, full_output=1)
    if convergence == 1:
        return sol
    return np.array([np.nan]*len(init))

def find_unique_equilibria(flow, starting_points):
    '''Return the list of unique equilibria of a flow 
    starting around starting_points'''
    equilibria = [] 
    roots = [findroot(flow, init) 
             for init in starting_points]
    # Only keep unique equilibria 
    for r in roots:
        if (not any(np.isnan(r)) and
            not any([all(np.isclose(r, x)) for x in equilibria])):
            equilibria.append(r)
    return equilibria

scenarios = list(np.linspace(0, 2.0, 5))
initial_conditions = [
    np.concatenate([[0]*192])
]
time = np.linspace(0, 400)

trajectory = {}
# trajectory = np.load('trajectories.npy', allow_pickle=True)
for i,param in enumerate(scenarios):
    print(param)
    for j,ic in enumerate(initial_conditions):
        # model.run(I0=1, R0=param, phases=EpidemicPhases(g0=1, itv0=0), verbose=False, attack=None, testing_parameters = TESTING_PARAMETERS_FROM_EPIDEMIOLOGY)
        # df = pd.read_csv('../output/' + 'result_data_BR.csv', sep=',')
        # trajectory[i,j] = df.values[:, 12:-18] # evolução do compartimento no tempo com R0 i e condição inicial j. T linhhas x VAR colunas
        func = Func(param)
        trajectory[i,j] = scipy.integrate.odeint(func.get_f(),
                                                 y0=ic,
                                                 t=time)

equilibria = {}
x = []
y = []
for i, param in enumerate(scenarios):
    # Find the position of the equilibirum around the endpoint of each trajectory. 
    func = Func(param)
    flow = func.get_f()
    starting_points = [trajectory[i,j][-1,:] for j 
                       in range(len(initial_conditions))] 
    equilibria = find_unique_equilibria(flow, starting_points)
    for eq in equilibria:
        x.append(param)
        I = []
        for count in range(age_strata):
            I.append(eq[count + 2])
        y.append(np.sum(I))
        print(eq)
    print('{} Equilibrium point(s) for parameters: {}'.format(len(equilibria), param))

plt.scatter(x, y)
plt.savefig('test2.png', format='png')