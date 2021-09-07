from .utils import read_matrix
import pandas as pd
import numpy as np
import csv

demographic_file = 'demographic_data.csv'
epidemiologic_file = 'epidemiology_data.csv'
testing_parameters_file = 'testing_parameters.csv'
hosp_uti_ratio = 0.3333

TESTING_PARAMETERS_FROM_EPIDEMIOLOGY = 1
TESTING_PARAMETERS_FROM_OWN_FILE = 2
TESTING_PARAMETERS_DIRECT = 3

class Parameters:

    def __init__(self, age_strata, num_phases, testing_parameters_mode, xI=None, xA=None, model=3, fatality=1.0, changeDir=False, scenario=''):
        if changeDir:
            self.input_folder = '../../input/cenarios/cenario' + scenario + '/'
        else:
            self.input_folder = ''

        self.num_phases = num_phases
        self.age_strata = age_strata
        self.model = model
        self.fatality = fatality

        self.data_epid = pd.read_csv(self.input_folder+epidemiologic_file).values[:, 1:].astype(np.float64)

        self.load_demographic_parameters()
        self.xI, self.xA = self.get_testing_parameters(testing_parameters_mode, xI, xA)
        self.load_epidemiology_parameters()
    
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
        self.gamma_QI = np.divide(np.multiply(self.xI, self.gamma_H + self.gamma_RI), (1 - self.xI))
        self.gamma_QA = np.divide(np.multiply(self.xA, self.gamma_RA), (1 - self.xA))
        self.tlc = np.divide(self.TC, self.phi)

        if self.model == 3:
            self.mu_cov = np.multiply(np.divide(self.fatality*self.TC, self.phi), np.multiply(1 - self.phi, self.gamma_H + self.gamma_RI))
        else:
            self.mu_cov = np.divide(np.multiply(gamma, self.fatality * self.TC), 1 - self.fatality * self.TC)
    
    def get_testing_parameters(self, testing_parameters_mode, input_xi, input_xa):
        age_strata = self.age_strata

        xI = np.zeros([self.num_phases, age_strata], dtype=np.float64)
        xA = np.zeros([self.num_phases, age_strata], dtype=np.float64)
        if testing_parameters_mode == TESTING_PARAMETERS_DIRECT:
            xI[0] = np.zeros(age_strata, dtype=np.float64)
            xA[0] = np.zeros(age_strata, dtype=np.float64)
            for k in range(1, self.num_phases):
                xI[k] = input_xi
                xA[k] = input_xa
        elif testing_parameters_mode == TESTING_PARAMETERS_FROM_OWN_FILE:
            data_testing_parameters = np.zeros([age_strata, age_strata], dtype=np.float64)
            n = 0
            with open(testing_parameters_file, "r") as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                next(spamreader, None)
                j = 0
                for row in spamreader:
                    for i in range(age_strata):
                        data_testing_parameters[j, i] = row[i+1]
                    j = j + 1
                n = int(j/2)
            for k in range(0, n):
                xI[k] = data_testing_parameters[2*k, :]
                xA[k] = data_testing_parameters[2*k+1, :]
            for k in range(n, self.num_phases):
                xI[k] = data_testing_parameters[2*(n-1), :]
                xA[k] = data_testing_parameters[2*(n-1)+1, :]
        else:  
            for k in range(0, self.num_phases):
                xI[k] = self.data_epid[2, :]
                xA[k] = self.data_epid[3, :]
        return xI, xA

    def load_demographic_parameters(self):
        data_demog = pd.read_csv(self.input_folder+demographic_file).values[:, 1:]
        self.pop = data_demog[0, :]
        self.mortalidade = data_demog[1, :] / (365 * 1000)
        self.natalidade = data_demog[2, :]

    def save(self, file):
        param_header = ['VARIAVEL', 'FAIXA_1', 'FAIXA_2', 'FAIXA_3', 'FAIXA_4', 'FAIXA_5', 'FAIXA_6', 'FAIXA_7', 'FAIXA_8',
                    'FAIXA_9', 'FAIXA_10', 'FAIXA_11', 'FAIXA_12', 'FAIXA_13', 'FAIXA_14', 'FAIXA_15', 'FAIXA_16']
        with open(file, 'w') as csvfile:
            spamwriter = csv.writer(csvfile)
            spamwriter.writerow(param_header)
            spamwriter.writerow(np.concatenate((['LAMBDA'], self.natalidade)))
            spamwriter.writerow(np.concatenate((['MORT_EQ'], self.mortalidade)))
            spamwriter.writerow(np.concatenate((['MORT_COV'], self.mu_cov)))
            spamwriter.writerow(np.concatenate((['THETA'], self.theta)))
            spamwriter.writerow(np.concatenate((['ALPHA'], self.alpha)))
            spamwriter.writerow(np.concatenate((['KSI'], self.ksi)))
            spamwriter.writerow(np.concatenate((['RHO'], self.rho)))
            spamwriter.writerow(np.concatenate((['PHI'], self.phi)))
            spamwriter.writerow(np.concatenate((['ETA'], self.eta)))
            spamwriter.writerow(np.concatenate((['A'], self.a)))
            spamwriter.writerow(np.concatenate((['GAMA_H'], self.gamma_H)))
            spamwriter.writerow(np.concatenate((['GAMA_HR'], self.gamma_HR)))
            spamwriter.writerow(np.concatenate((['GAMA_RI'], self.gamma_RI)))
            spamwriter.writerow(np.concatenate((['GAMA_RA'], self.gamma_RA)))
            spamwriter.writerow(np.concatenate((['GAMA_RQI'], self.gamma_RQI)))
            spamwriter.writerow(np.concatenate((['GAMA_RQA'], self.gamma_RQA)))
            spamwriter.writerow(np.concatenate((['GAMA_HQI'], self.gamma_HQI)))
            spamwriter.writerow(np.concatenate((['TC'], self.TC)))
            spamwriter.writerow(np.concatenate((['TLC'], self.tlc)))

    def read_gamma_H(self):
        return self.gamma_H