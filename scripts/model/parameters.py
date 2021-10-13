import pandas as pd
import numpy as np
import csv
import cmodels
from .contact_matrices import ContactMatrixFactory
import sys
from scipy import linalg

demographic_file = 'demographic_data.csv'
epidemiologic_file = 'epidemiology_data.csv'
testing_parameters_file = 'testing_parameters.csv'
hosp_uti_ratio = 0.3333

E = 0
A = 1
I = 2
Qi = 3
Qa = 4
H = 5

class EpidemicPhases:
    def __init__(self, g0, itv0, xI0, xA0):
        self.g = [g0]
        self.days = [0]
        self.itv_level = [itv0]
        self.xIs = [xI0]
        self.xAs = [xA0]

    def add_phase(self, g, start_day, itv_level, xI, xA):
        self.g.append(g)
        self.days.append(start_day)
        self.itv_level.append(itv_level)
        self.xIs.append(xI)
        self.xAs.append(xA)

    def get_length(self):
        return len(self.g)

class Parameters:

    def __init__(self, scenario, age_strata, phases, model=3, fatality=1.0, I0=1):
        self.input_folder = '../input/cenarios/cenario' + scenario + '/'

        self.I0 = I0
        self.scenario = scenario
        self.age_strata = age_strata
        self.model = model
        self.fatality = fatality

        self.data_epid = pd.read_csv(self.input_folder+epidemiologic_file).values[:, 1:].astype(np.float64)

        self.pop, self.mu_eq, self.natalidade, self.rel_pop = self.load_demographic_parameters()
        self.load_epidemiology_parameters()
        self.cparameters = self.load_cparameters()

        self.R0, self.dynamic_parameters = self.load_dynamic_parameters(phases, self.rel_pop, fatality, age_strata)

        self.get_initial_matrix(self.R0, I0)
        
    def load_epidemiology_parameters(self):
        age_strata = self.age_strata
        
        tlc = np.zeros(age_strata, dtype=np.float64)
        
        # Read epidemiologic data
        self.theta = np.zeros(age_strata, dtype=np.float64)
        self.TC = self.data_epid[0, :]
        self.phi = self.data_epid[1, :]
        dI = self.data_epid[4, :]
        dL = self.data_epid[5, :]
        dH = self.data_epid[6, :]
        dA = self.data_epid[7, :]
        self.eta = self.data_epid[8, :]
        self.rho = self.data_epid[9, :]
        self.alpha = self.data_epid[10, :]
        self.ksi = self.data_epid[11, :]

        self.eta = hosp_uti_ratio * self.phi
        self.a = 1 - np.exp(-np.divide(1, dL))
        self.gamma = 1 - np.exp(-np.divide(1, dI))
        self.gamma_RI = 1 - np.exp(-np.divide(1, dI))
        self.gamma_RA = 1 - np.exp(-np.divide(1, dI))
        self.gamma_RQI = 1 - np.exp(-np.divide(1, dI))
        self.gamma_RQA = 1 - np.exp(-np.divide(1, dI))
        self.gamma_HR = 1 - np.exp(-np.divide(1, dA))
        self.gamma_H = np.divide(np.multiply(self.gamma_RI, self.phi), 1 - self.phi)
        self.gamma_HQI = np.divide(np.multiply(self.gamma_RQI, self.phi), 1 - self.phi)
        self.tlc = np.divide(self.TC, self.phi)

        self.mu_cov = self.load_mu_cov()
    
    def load_mu_cov(self):
        if self.model == 3:
            return np.multiply(np.divide(self.fatality*self.TC, self.phi), np.multiply(1 - self.phi, self.gamma_H + self.gamma_RI))
        else:
            return np.divide(np.multiply(gamma, self.fatality * self.TC), 1 - self.fatality * self.TC)

    def load_demographic_parameters(self):
        data_demog = pd.read_csv(self.input_folder+demographic_file).values[:, 1:]
        pop = data_demog[0, :]
        mu_eq = data_demog[1, :] / (365 * 1000)
        natalidade = data_demog[2, :]
        rel_pop = np.divide(pop, np.sum(pop))
        return pop, mu_eq, natalidade, rel_pop

        
    def get_initial_matrix(self, R0, I0):
        self.initial = np.zeros((13, 16))
        I_0_vec = np.zeros(self.pop.size)
        I_0_vec[13] = I0
        self.initial[0, :] = self.pop
        self.initial[1, :] = I0 * R0 * self.rel_pop
        self.initial[2, :] = I_0_vec
        self.initial[4, :] = I_0_vec
        self.initial[5, :] = np.zeros((16, ))
        self.initial[6, :] = np.zeros((16, ))
        self.initial[9, :] = self.pop
    
    def load_cparameters(self):
        cparameters = cmodels.StaticParameters()
        cparameters.Lambda = self.natalidade
        cparameters.mu_eq = self.mu_eq
        cparameters.alpha = self.alpha
        cparameters.ksi = self.ksi
        cparameters.rho = self.rho
        cparameters.phi = self.phi
        cparameters.eta = self.eta
        cparameters.a = self.a
        cparameters.mu_cov = self.mu_cov
        cparameters.theta = self.theta
        cparameters.gama_A = np.zeros((16, 1))
        cparameters.gama_H = self.gamma_H
        cparameters.gama_HR = self.gamma_HR
        cparameters.gama_RI = self.gamma_RI
        cparameters.gama_RA = self.gamma_RA
        cparameters.gama_RQI = self.gamma_RQI
        cparameters.gama_RQA = self.gamma_RQA
        cparameters.gama_HQI = self.gamma_HQI
        cparameters.Tc = self.TC
        cparameters.gama = self.gamma
        cparameters.Tlc = self.tlc
        cparameters.rel_pop = self.rel_pop
        return cparameters

    def update_phases(self, phases, attack=1.0):
        self.fatality = attack
        # updates mortality
        self.cparameters.mu_cov = self.load_mu_cov()
        self.R0, self.dynamic_parameters = self.load_dynamic_parameters(phases, self.rel_pop, self.fatality, self.age_strata)
        self.get_initial_matrix(self.R0, self.I0)
        return self.R0

    def load_dynamic_parameters(self, phases, rel_pop, attack, age_strata):
        prob_t = phases.g
        num_phases = phases.get_length()
        itv_list = phases.itv_level
        days = phases.days

        cmf = ContactMatrixFactory(age_strata, self.scenario)

        contact_matrix = cmf.get_matrix(itv_list[0], rel_pop, attack, age_strata)
        dynamic_parameters = cmodels.DynamicParameters(phases.xIs[0], phases.xAs[0], prob_t[0] * contact_matrix)
        R0 = cmodels.calculateR0(self.cparameters, dynamic_parameters)
        # Other phases
        for i in range(1, num_phases):
            contact_matrix = cmf.get_matrix(itv_list[i], rel_pop, attack, age_strata)
            dynamic_parameters.addPhase(
                days[i], 
                phases.xIs[i],
                phases.xAs[i],
                prob_t[i] * contact_matrix
            )
        return R0, dynamic_parameters

    def calculateR0(self, xI, xA, betas):
        age_strata = self.age_strata
        self.gamma_QI = np.multiply(xI, self.gamma_H+self.gamma_RI) / (1-xI)
        self.gamma_QA = np.multiply(xA, self.gamma_RA) / (1-xA)

        # A unica entrada de novas infecções é em E, vindo de A, I ou E
        F = np.zeros((6*age_strata, 6*age_strata))

        # Some matrix manipulation
        B_1 = np.multiply(self.ksi, np.multiply(betas, np.tile(self.rel_pop.reshape(age_strata, 1), (1, age_strata))))
        B_2 = np.multiply(self.alpha, np.multiply(betas, np.tile(self.rel_pop.reshape(age_strata, 1), (1, age_strata))))
        B_3 = np.multiply(betas, np.tile(self.rel_pop.reshape(age_strata, 1), (1, age_strata)))
        F[E * age_strata: (E + 1) * age_strata, E * age_strata: (E + 1) * age_strata] = B_1
        F[E * age_strata: (E + 1) * age_strata, A * age_strata: (A + 1) * age_strata] = B_2
        F[E * age_strata: (E + 1) * age_strata, I * age_strata: (I + 1) * age_strata] = B_3

        # Transições entre compartimentos de infectados
        V = np.zeros((6*age_strata, 6*age_strata))
        # Exposto -> Exposto
        for i in range(age_strata):
            j = E + i
            V[j, j] = self.mu_eq[i] + self.a[i]
        # Assintomático
        for i in range(age_strata):
            j = A*age_strata + i
            V[j, j] = self.mu_eq[i] + self.gamma_RA[i] + self.gamma_QA[i]
            V[j, E*age_strata+i] = -self.a[i] * (1-self.rho[i])
        # Infectado sintomático
        for i in range(age_strata):
            j = I*age_strata + i
            V[j, j] = self.mu_eq[i] + self.theta[i]*self.mu_cov[i] + self.gamma_H[i] + self.gamma_RI[i] + self.gamma_QI[i]
            V[j, E*age_strata+i] = -self.a[i]*self.rho[i]
        # Quarentenado sintomático
        for i in range(age_strata):
            j = Qi*age_strata + i
            V[j, j] = self.mu_eq[i] + self.gamma_HQI[i] + self.gamma_RQI[i]
            V[j, I*age_strata+i] = -self.gamma_QI[i]
        # Quarentenado assintomático
        for i in range(age_strata):
            j = Qa*age_strata + i
            V[j, j] = self.mu_eq[i] + self.gamma_RQA[i]
            V[j, A*age_strata+i] = -self.gamma_QA[i]
        # Hospitalizado
        for i in range(age_strata):
            j = H*age_strata + i
            V[j, j] = self.mu_eq[i] + self.gamma_HR[i] + self.mu_cov[i]
            V[j, I*age_strata+i] = -self.gamma_H[i]
            V[j, Qi*age_strata+i] = -self.gamma_HQI[i]

        # print(V)
        v_1 = np.linalg.inv(V)
        np.set_printoptions(threshold=sys.maxsize)
        np.savetxt("V.csv", V, delimiter=",")
        np.savetxt("V-1.csv", v_1, delimiter=",")
        np.savetxt("F.csv", F, delimiter=",")
        next_gen_matrix = np.matmul(F, v_1)
        w, _ = linalg.eig(next_gen_matrix)
        R0 = np.max(np.abs(w))
        # if verbose:
        print("R0 obtido: "+str(R0))
        return R0
            
