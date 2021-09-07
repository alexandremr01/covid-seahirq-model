# coding: utf-8
import numpy as np
import csv
import os
import sys
import pandas as pd
from .nextgen import calculateR0
from .utils import read_matrix
from .parameters import Parameters
from .contact_matrices import ContactMatrixFactory

TESTING_PARAMETERS_FROM_EPIDEMIOLOGY = 1
TESTING_PARAMETERS_FROM_OWN_FILE = 2
TESTING_PARAMETERS_DIRECT = 3

def scen_gen(input_folder, model, I_0, R0, phases, fatality=1.0, verbose=False, input_xi=None, input_xa=None, attack = None,
                testing_parameters_mode=TESTING_PARAMETERS_FROM_EPIDEMIOLOGY):
    
    age_strata = 16
    prev_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join('..', 'input', 'cenarios', input_folder))

    prob_t = phases.g
    num_phases = phases.get_length()
    itv_list = phases.itv_level
    days = phases.days

    ## show values ##
    if verbose:
        print("Input file: %s" % input_folder)
        print("Days: %s" % days)
        print("Model: %s" % model)
        print("I0:%s" % I_0)
        print("R0: %s" % R0)
        print("Correction post: %s" % prob_t)
        print("Intervention cases: %s" % itv_list)
        print("N: %s" % num_phases)

        print(u'Diretório inicial ' + prev_dir)
        print(u'Movendo para o diretório de entrada (input) ')
        print(u'Diretório de entrada (input) ' + os.path.abspath(os.curdir))

    # Read and save parameters
    parameters = Parameters(age_strata, num_phases, testing_parameters_mode, input_xi, input_xa)
    parameters.save('parameters.csv')

    pop = parameters.pop
    rel_pop = np.divide(pop, np.sum(pop))  # partitioning of the cases, for I0 > age_strata

    # Choose intervention matrixes
    contact_matrices = []
    contactMatrixFactory = ContactMatrixFactory(age_strata)

    # Phase 0: get "actual" R0
    contact_matrix = contactMatrixFactory.get_matrix(itv_list[0], rel_pop, attack, age_strata, input_folder, parameters, R0)
    contact_matrices.append( prob_t[0] * contact_matrix )
    actual_R0 = calculateR0(input_folder, contact_matrix, parameters.xI[0], parameters.xA[0], False)

    # Other phases
    for i in range(1, num_phases):
        contact_matrix = contactMatrixFactory.get_matrix(itv_list[i], rel_pop, attack, age_strata, input_folder, parameters, R0)
        contact_matrices.append( prob_t[i] * contact_matrix )
    
    write_initial(pop, actual_R0, I_0, rel_pop)
    write_beta_gama(num_phases, age_strata, days, parameters, contact_matrices)
    os.chdir(os.path.join('..', '..', '..', 'scripts'))
    return actual_R0

def write_initial(pop, R0, I_0, rel_pop):
    I_0_vec = np.zeros(pop.size)
    I_0_vec[13] = I_0

    with open('initial.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(np.concatenate((['POPULATION_S'], pop)))
        spamwriter.writerow(np.concatenate((['EXPOSED_E'], I_0 * R0 * rel_pop)))
        spamwriter.writerow(np.concatenate((['ASYMPTOMATIC_A'], I_0_vec)))
        spamwriter.writerow(np.concatenate((['INFECTED_I'], I_0_vec)))
        spamwriter.writerow(np.concatenate((['RECOVERED_R'], np.zeros(pop.size))))
        spamwriter.writerow(np.concatenate((['RECOVERED_SYMPTOMATIC_Ri'], np.zeros(pop.size))))

def write_beta_gama(num_phases, age_strata, days, parameters, matrix):
    beta_gama_header = ['DAY', 'GAMA_F1', 'GAMA_F2', 'GAMA_F3', 'GAMA_F4', 'GAMA_F5', 'GAMA_F6', 'GAMA_F7', 'GAMA_F8',
                        'GAMA_F9', 'GAMA_F10', 'GAMA_F11', 'GAMA_F12', 'GAMA_F13', 'GAMA_F14', 'GAMA_F15', 'GAMA_F16',
                        'xI_F1', 'xI_F2', 'xI_F3', 'xI_F4', 'xI_F5', 'xI_F6', 'xI_F7',
                        'xI_F8', 'xI_F9', 'xI_F10', 'xI_F11', 'xI_F12', 'xI_F13',
                        'xI_F14', 'xI_F15', 'xI_F16',
                        'xA_F1', 'xA_F2', 'xA_F3', 'xA_F4', 'xA_F5', 'xA_F6', 'xA_F-',
                        'xA_F8', 'xA_F9', 'xA_F10', 'xA_F11', 'xA_F12', 'xA_F13',
                        'xA_F14', 'xA_F15', 'xA_F16',
                        'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX',
                        'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX',
                        'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX', 'BETA_MATRIX']

    with open('beta_gama.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(beta_gama_header)
        for k in range(num_phases):
            for i in range(0, age_strata):
                if i == 0:
                    spamwriter.writerow(np.concatenate(([days[k]], parameters.gamma, parameters.xI[k], parameters.xA[k], matrix[k][i, :])))
                else:
                    spamwriter.writerow(np.concatenate(([days[k]], ['']*48, matrix[k][i, :])))
            spamwriter.writerow(['']*65)

