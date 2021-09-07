# coding: utf-8
import numpy as np
import csv
import os
import sys
import argparse
from scipy.sparse.linalg import eigs

nu_exp = 0.4  # exposed contamination
hosp_uti_ratio = 0.3333


#### Some function definition
def symmetrize(matrix_C, PP, age_n):
    c_m = np.zeros([age_n, age_n])
    for ii in range(0, age_n):
        for jj in range(0, age_n):
            c_m[ii, jj] = 0.5 * (matrix_C[ii, jj] * PP[ii] / PP[jj] + matrix_C[jj, ii] * PP[jj] / PP[ii])
    return c_m


#### Pick the chosen scenario matrix
def pick_intervention(itv_ic, *arg):
    if itv_ic == 0:
        return arg[0]
    elif itv_ic == 1:
        return arg[1]
    elif itv_ic == 2:
        return arg[2]
    elif itv_ic == 3:
        return arg[3]
    elif itv_ic == 4:
        return arg[4]
    elif itv_ic == 5:
        return arg[5]
    elif itv_ic == 6:
        return arg[6]
    elif itv_ic == 7:
        return arg[7]
    elif itv_ic == 8:
        return arg[8]
    elif itv_ic == 9:
        return arg[9]
    elif itv_ic == 10:
        return arg[10]
    elif itv_ic == 11:
        return arg[11]
    else:
        return arg[0]


##### Process command line options
##### Variable parameters, for error estimation within reasonable bounds
parser = argparse.ArgumentParser(description='This script computes the scenario matrices and other inputs for '
                                             'csv_to_input.')
parser.add_argument('-i', '--input_file', help='Input file name', required=True)
parser.add_argument('-d', '--day', type=float, nargs='+', help='Days of measure beginning - four values required',
                    required=True)
parser.add_argument('-m', '--model', type=int, help='(0) SIR, (1) SEIR, (2) SEAIR, (3) SEAHIR - default is SEAHIR ',
                    required=True)
parser.add_argument('-I0', '--I0_value', type=float, help='Number of initial infected', required=True)
parser.add_argument('-R0', '--R0_value', type=float, help='Reproduction number', required=True)
# parser.add_argument('-Rp', '--R0_post', type=float, help='Post outbreak reproduction number', required=False)
# parser.add_argument('-n', '--num_stages', type=int, help='Number of stages of the desease evolution', required=True)
parser.add_argument('-z', '--zeta', type=float, help='Ratio between deaths of infecteds that were never hospitalized and of those who were.', required=True)
parser.add_argument('-t', '--tau', type=float, help='Ratio between phi and IFR.', required=True)
parser.add_argument('-n', '--num_phases', type=int, help='Number of phases of interventions/testing.', required=False, default=4)
parser.add_argument('-p', '--prob', type=float, nargs='+', help='Overall transmission probability decrease',
                    required=True) 
parser.add_argument('-itv', '--intervention', type=int, nargs='+', help='Four integer numbers, each representing an '
                                                                      'intervention case. See README file for '
                                                                      'intervention table.',
                    required=True)
parser.add_argument('-tp', '--test_parameters', action='store_true', help='If set, this argument reads xI and xA from'
                                                                         'test_parameters.csv instead of epidemiology_data.csv',
                    required=False)
parser.add_argument('-f', '--fatality', type=float, help='Factor that multiplies uniformly the fatality '
                                                         'rate. Use to estimate uncertainty. ',
                    required=False)
parser.add_argument('-s', action='store_true', help='If set, this argument print the equivalent Rt of each scenario',
                    required=False)
parser.add_argument('-ex', type=int, nargs=2, help='Extra intervention case. Require the intervention code and '
                                                   'then the date. See README file for intervention table.',
                    required=False)
parser.add_argument('-ifrs', '--ifrs', type=float, nargs=16, help='Age specific ifr.', required=False)

args = parser.parse_args()

## show values ##
print ("Input file: %s" % args.input_file)
print ("Days: %s" % args.day)
print ("Model: %s" % args.model)
print ("I0:%s" % args.I0_value)
print ("R0: %s" % args.R0_value)
print ("Correction post: %s" % args.prob)
print ("Intervention cases: %s" % args.intervention)
print ("Tau: %s" % args.tau)
print ("Zeta: %s" % args.zeta)
print ("N: %s" % args.num_phases)

if args.fatality is None:
    f_fac = 1.0
else:
    f_fac = float(args.fatality)

I_0 = int(args.I0_value)
R0 = float(args.R0_value)
R0_post = R0 #float(args.R0_post)
model = int(args.model)
days = [float(x) for x in args.day]
tau = float(args.tau)
zeta = float(args.zeta)
itv_list = [int(x) for x in args.intervention]
input_folder = args.input_file
prob_t = np.array([float(x) for x in args.prob], dtype=np.float64)# protection decrease of transmission probability

num_phases = args.num_phases
age_strata = 16
t_days = 400
demographic_file = 'demographic_data.csv'
epidemiologic_file = 'epidemiology_data.csv'
contact_matrix_all_file = 'contact_matrix_all.csv'
contact_matrix_home_file = 'contact_matrix_home.csv'
contact_matrix_work_file = 'contact_matrix_work.csv'
contact_matrix_school_file = 'contact_matrix_school.csv'
contact_matrix_other_file = 'contact_matrix_other.csv'
testing_parameters_file = 'testing_parameters.csv'

matrix_all = np.zeros([age_strata, age_strata], dtype=np.float64)
matrix_home = np.zeros([age_strata, age_strata], dtype=np.float64)
matrix_work = np.zeros([age_strata, age_strata], dtype=np.float64)
matrix_school = np.zeros([age_strata, age_strata], dtype=np.float64)
matrix_other = np.zeros([age_strata, age_strata], dtype=np.float64)

demog_variables = 3
data_demog = np.zeros([demog_variables, age_strata], dtype=np.float64)
pop = np.zeros(age_strata, dtype=np.float64)
mortalidade = np.zeros(age_strata, dtype=np.float64)
natalidade = np.zeros(age_strata, dtype=np.float64)

epid_variables = 12
data_epid = np.zeros([epid_variables, age_strata], dtype=np.float64)
TC = np.zeros(age_strata, dtype=np.float64)
phi = np.zeros(age_strata, dtype=np.float64)
dI = np.zeros(age_strata, dtype=np.float64)
dL = np.zeros(age_strata, dtype=np.float64)
dH = np.zeros(age_strata, dtype=np.float64)
dA = np.zeros(age_strata, dtype=np.float64)
eta = np.zeros(age_strata, dtype=np.float64)
pho = np.zeros(age_strata, dtype=np.float64)
alpha = np.zeros(age_strata, dtype=np.float64)
tlc = np.zeros(age_strata, dtype=np.float64)
theta = np.zeros(age_strata, dtype=np.float64)

# Testagem
xI = np.zeros([num_phases, age_strata], dtype=np.float64)
xA = np.zeros([num_phases, age_strata], dtype=np.float64)


curr_dir = os.path.abspath(os.curdir)
print(u'Diretório atual ' + curr_dir)
print(u'Movendo para o diretório de entrada (input) ')
os.chdir("..")
os.chdir(os.path.join('input', 'cenarios', input_folder))

curr_dir = os.path.abspath(os.curdir)
print(u'Diretório de entrada (input) ' + curr_dir)

# Read demographic data
with open(demographic_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            data_demog[j, i] = row[i + 1]
        j = j + 1

pop = data_demog[0, :]
mortalidade = data_demog[1, :] / (365 * 1000)
natalidade = data_demog[2, :]

# Read epidemiologic data
with open(epidemiologic_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            data_epid[j, i] = row[i + 1]
        j = j + 1

# Read contact matrices
# All
with open(contact_matrix_all_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            matrix_all[j, i] = row[i]
        j = j + 1

# Home
with open(contact_matrix_home_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            matrix_home[j, i] = row[i]
        j = j + 1

# Work
with open(contact_matrix_work_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            matrix_work[j, i] = row[i]
        j = j + 1

# School
with open(contact_matrix_school_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            matrix_school[j, i] = row[i]
        j = j + 1

# Other
with open(contact_matrix_other_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            matrix_other[j, i] = row[i]
        j = j + 1

if args.ifrs is None:
    TC = data_epid[0, :]
    phi = tau*TC
else:
    TC = np.array([float(x) for x in args.ifrs], dtype=np.float64)
    phi = tau*TC

#theta = zeta * phi
theta = zeta * np.ones(age_strata)
print("THETA: ")
print(theta)
dI = data_epid[4, :]
dL = data_epid[5, :]
dH = data_epid[6, :]
dA = data_epid[7, :]
eta = data_epid[8, :]
rho = data_epid[9, :]
alpha = data_epid[10, :]
#zeta = data_epid[11, :]

# Testagem
if args.test_parameters: #Pega de testing_parameters
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
    print("n="+str(n))
    print(data_testing_parameters[4])
    for k in range(0, n):
        xI[k] = data_testing_parameters[2*k, :]
        xA[k] = data_testing_parameters[2*k+1, :]
    for k in range(n, num_phases):
        xI[k] = data_testing_parameters[2*(n-1), :]
        xA[k] = data_testing_parameters[2*(n-1)+1, :]
else:  # Pega de epidemiology_data
    print("Pegando de epidemiology")
    for k in range(0, num_phases):
        xI[k] = data_epid[2, :]
        xA[k] = data_epid[3, :]


eta = hosp_uti_ratio * phi

# calcula os valores de 'a' e gama's e mu_cov

a = 1 - np.exp(-np.divide(1, dL))
gamma = 1 - np.exp(-np.divide(1, dI))
gamma_RI = 1 - np.exp(-np.divide(1, dI))
gamma_RA = 1 - np.exp(-np.divide(1, dI))
gamma_RQI = 1 - np.exp(-np.divide(1, dI))
gamma_RQA = 1 - np.exp(-np.divide(1, dI))
gamma_HR = 1 - np.exp(-np.divide(1, dA))
gamma_H = np.divide(np.multiply(gamma_RI, phi), 1 - phi)
gamma_HQI = np.divide(np.multiply(gamma_RQI, phi), 1 - phi)
gamma_QI = np.divide(np.multiply(xI, gamma_H + gamma_RI), (1 - xI))
gamma_QA = np.divide(np.multiply(xA, gamma_RA), (1 - xA))
tlc = np.divide(TC, phi)

if model == 3:
    #mu_cov = np.multiply(np.divide(f_fac*TC, phi), gamma_H + gamma_RI)
    mu_cov_num = np.multiply(TC, gamma_RI)
    mu_cov_den = np.multiply(theta+phi, 1-phi)
    mu_cov = np.divide(mu_cov_num, mu_cov_den)
else:
    mu_cov = np.divide(np.multiply(gamma, f_fac * TC), 1 - f_fac * TC)

# cria arquivo initial.csv

rel_pop = np.divide(pop, np.sum(pop))  # partitioning of the cases, for I0 > age_strata
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
# cria arquivo parameters.csv

param_header = ['VARIAVEL', 'FAIXA_1', 'FAIXA_2', 'FAIXA_3', 'FAIXA_4', 'FAIXA_5', 'FAIXA_6', 'FAIXA_7', 'FAIXA_8',
                'FAIXA_9', 'FAIXA_10', 'FAIXA_11', 'FAIXA_12', 'FAIXA_13', 'FAIXA_14', 'FAIXA_15', 'FAIXA_16']

with open('parameters.csv', 'w') as csvfile:
    spamwriter = csv.writer(csvfile)
    spamwriter.writerow(param_header)
    spamwriter.writerow(np.concatenate((['LAMBDA'], natalidade)))
    spamwriter.writerow(np.concatenate((['MORT_EQ'], mortalidade)))
    spamwriter.writerow(np.concatenate((['MORT_COV'], mu_cov)))
    spamwriter.writerow(np.concatenate((['THETA'], theta)))
    spamwriter.writerow(np.concatenate((['ALPHA'], alpha)))
    spamwriter.writerow(np.concatenate((['RHO'], rho)))
    spamwriter.writerow(np.concatenate((['PHI'], phi)))
    spamwriter.writerow(np.concatenate((['ETA'], eta)))
    spamwriter.writerow(np.concatenate((['A'], a)))
    spamwriter.writerow(np.concatenate((['GAMA_H'], gamma_H)))
    spamwriter.writerow(np.concatenate((['GAMA_HR'], gamma_HR)))
    spamwriter.writerow(np.concatenate((['GAMA_RI'], gamma_RI)))
    spamwriter.writerow(np.concatenate((['GAMA_RA'], gamma_RA)))
    spamwriter.writerow(np.concatenate((['GAMA_RQI'], gamma_RQI)))
    spamwriter.writerow(np.concatenate((['GAMA_RQA'], gamma_RQA)))
    spamwriter.writerow(np.concatenate((['GAMA_HQI'], gamma_HQI)))
    # spamwriter.writerow(np.concatenate((['GAMA_QI'], gamma_QI)))
    # spamwriter.writerow(np.concatenate((['GAMA_QA'], gamma_QA)))
    spamwriter.writerow(np.concatenate((['TC'], TC)))
    spamwriter.writerow(np.concatenate((['TLC'], tlc)))

# cria arquivos beta_gama.csv

# Escala a matriz de contato pelo perfil demográfico e escalona a mesma pelo auto-valor dominante

P = np.outer(pop, np.divide(1, pop))
# C_sym_all = 0.5*(np.dot(np.transpose(P), matrix_all) + np.dot(np.transpose(matrix_all), P))
C_sym_home = 0.5 * (matrix_home + symmetrize(matrix_home, pop, age_strata))
C_sym_school = 0.5 * (matrix_school + symmetrize(matrix_school, pop, age_strata))
C_sym_work = 0.5 * (matrix_work + symmetrize(matrix_work, pop, age_strata))
C_sym_other = 0.5 * (matrix_other + symmetrize(matrix_other, pop, age_strata))
C_sym = C_sym_home + C_sym_work + C_sym_school + C_sym_other
for i in range(0, age_strata):
    for j in range(0, age_strata):
        C_sym[i, j] = C_sym[i, j] * pop[i] / pop[j]
        C_sym_home[i, j] = C_sym_home[i, j] * pop[i] / pop[j]
        C_sym_school[i, j] = C_sym_school[i, j] * pop[i] / pop[j]
        C_sym_work[i, j] = C_sym_work[i, j] * pop[i] / pop[j]
        C_sym_other[i, j] = C_sym_other[i, j] * pop[i] / pop[j]
w, v = eigs(C_sym)
# Main eigenvector is normalized, norm_vec = 1. Mean square_vec=1/age_strata
# print(w.max())
# eig_vec = v[:, 0]
# print(v[:, 0])
# norm_vec = np.dot(eig_vec, eig_vec)
# square_vec = np.multiply(eig_vec, eig_vec)
# print(np.mean(norm_vec), 1./16.)
eig_value = np.real(w.max())
if model == 2 or 3:
    beta = R0 * gamma[0] / (rho.max() + np.mean(alpha) * (1 - rho.max())) / eig_value
    beta_val = R0 * gamma[0] / (rho.max() + np.mean(alpha) * (1 - rho.max()))
elif model == 3:
    beta = R0 * gamma[0] / (rho.max() + np.mean(alpha) * (1 - rho.max()) + gamma[0] * nu_exp / a[0]) / eig_value
else:
    beta = R0 * gamma[0]
C_all_pre = beta * C_sym * age_strata  ## itv_id = 0
C_home_pre = beta * C_sym_home * age_strata
C_work_pre = beta * C_sym_work * age_strata
C_school_pre = beta * C_sym_school * age_strata
C_other_pre = beta * C_sym_other * age_strata

if model == 2 or model == 3:
    beta = R0_post * gamma[0] / (rho.max() + np.mean(alpha) * (1 - rho.max())) / eig_value
    beta_val = R0_post * gamma[0] / (rho.max() + np.mean(alpha) * (1 - rho.max()))
elif model == 3:
    beta = R0 * gamma[0] / (rho.max() + np.mean(alpha) * (1 - rho.max()) + gamma[0] * nu_exp / a[0]) / eig_value
else:
    beta = R0 * gamma[0]

C_home_post = beta * C_sym_home * age_strata
C_work_post = beta * C_sym_work * age_strata
C_school_post = beta * C_sym_school * age_strata
C_other_post = beta * C_sym_other * age_strata
C_all_post = C_home_post + C_work_post + C_school_post + C_other_post  ## itv_id = 1

# Build matrix for scenarios
I_old = np.diag(np.ones(age_strata))
A_home = np.diag(np.ones(age_strata))
B_other = np.diag(np.ones(age_strata))
B_school = np.diag(np.ones(age_strata))
B_lock = np.diag(np.ones(age_strata))
B_strong_lock = np.diag(np.ones(age_strata))
W_work = np.diag(np.ones(age_strata))
W_lock = np.diag(np.ones(age_strata))
W_strong_lock = np.diag(np.ones(age_strata))
for i in range(age_strata - 5, age_strata):
    I_old[i, i] = 0.1
for i in range(0, 4):
    A_home[i, i] = 1.5
for i in range(4, age_strata):
    A_home[i, i] = 1.1
for i in range(0, 4):
    B_school[i, i] = 0.4
for i in range(0, 4):
    B_other[i, i] = 0.4
for i in range(4, age_strata):
    B_other[i, i] = 0.6
for i in range(0, 4):
    B_lock[i, i] = 0.3
for i in range(4, age_strata):
    B_lock[i, i] = 0.5
for i in range(0, 4):
    B_strong_lock[i, i] = 0.1
for i in range(4, age_strata):
    B_strong_lock[i, i] = 0.1
for i in range(0, age_strata):
    W_work[i, i] = 0.5
for i in range(0, age_strata):
    W_lock[i, i] = 0.4
for i in range(0, age_strata):
    W_strong_lock[i, i] = 0.3

# Fechamento de escolas apenas
C_all_school = np.dot(A_home, C_home_post) + C_work_post + np.dot(B_school, C_other_post)

# Fechamento de escolas e distanciamento social
C_all_school_other = np.dot(A_home, C_home_post) + C_work_post + + np.dot(B_other, C_other_post)

# Fechamento de escolas, distanciamento social e trabalho
C_all_school_other_work = np.dot(A_home, C_home_post) + np.dot(W_work, C_work_post) + np.dot(B_other, C_other_post)

# Lockdown sem isolamento de idoso
C_all_lock = np.dot(A_home, C_home_post) + np.dot(W_lock, C_work_post) + np.dot(B_lock, C_other_post)

C_work_old = np.dot(np.dot(I_old, C_work_post), I_old)
C_school_old = np.dot(np.dot(I_old, C_school_post), I_old)
C_other_old = np.dot(np.dot(I_old, C_other_post), I_old)
# Isolamento dos idosos em, com redução de X% dos seus contatos
C_all_old = C_home_post + C_work_old + C_school_old + C_other_old

# Isolamento dos idosos, fechamento de escolas
C_all_old_school = np.dot(A_home, C_home_post) + C_work_old + np.dot(B_school, C_other_old)

# Isolamento dos idosos, fechamento de escolas e distanciamento social
C_all_old_school_other = np.dot(A_home, C_home_post) + C_work_old + np.dot(B_other, C_other_old)

# Isolamento dos idosos, fechamento de escolas, distanciamento social e distanciamento no trabalho
C_all_old_school_other_work = np.dot(A_home, C_home_post) + np.dot(W_work, C_work_old) + np.dot(B_other, C_other_old)

# Lockdown com isolamento de idoso
C_all_old_lock = np.dot(A_home, C_home_post) + np.dot(W_lock, C_work_old) + np.dot(B_lock, C_other_old)

# Lockdown com isolamento de idoso
C_all_old_strong_lock = np.dot(A_home, C_home_post) + np.dot(W_strong_lock, C_work_old) + np.dot(B_strong_lock,
                                                                                                 C_other_old)

# Escolhe as matrizes de intervenção
matrix = []
for i in range(num_phases):
    matrix.append( prob_t[i] * pick_intervention(itv_list[i], C_all_pre, C_all_post, C_all_school, C_all_school_other,
                                         C_all_school_other_work, C_all_lock, C_all_old, C_all_old_school,
                                         C_all_old_school_other, C_all_old_school_other_work, C_all_old_lock,
                                         C_all_old_strong_lock) )

if args.ex is not None:
    opt_ex = int(args.ex[0])
    day_next_4 = int(args.ex[1])
    matrix_4 = prob_t[3] * pick_intervention(opt_ex, C_all_pre, C_all_post, C_all_school, C_all_school_other,
                                             C_all_school_other_work, C_all_lock, C_all_old, C_all_old_school,
                                             C_all_old_school_other, C_all_old_school_other_work, C_all_old_lock,
                                             C_all_old_strong_lock)

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
space_16 = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']

space_48 = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
            '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
            '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']

with open('beta_gama.csv', 'w') as csvfile:
    spamwriter = csv.writer(csvfile)
    spamwriter.writerow(beta_gama_header)
    for k in range(num_phases):
        for i in range(0, age_strata):
            if i == 0:
                spamwriter.writerow(np.concatenate(([days[k]], gamma, xI[k], xA[k], matrix[k][i, :])))
            else:
                spamwriter.writerow(np.concatenate(([days[k]], space_48, matrix[k][i, :])))
        spamwriter.writerow(np.concatenate(([''], space_48, space_16)))

print(u'Voltando para o diretório de script')
os.chdir("..")
os.chdir("..")
os.chdir("..")
os.chdir("scripts")

if args.s is True:
    R = []
    for k in range(num_phases):
        w, v = eigs(matrix[k])
        Rk = (np.real(w.max()) / age_strata) * (rho.max() + np.mean(alpha) * (1 - rho.max()) + gamma[0] * nu_exp / a[0]) / \
             gamma[0]
        R.append(Rk)
        print('R'+str(k)+' value is: ', Rk)
    os.chdir("..")
    os.chdir(os.path.join('input', 'cenarios', input_folder))
    with open('optimized_parameters.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(np.concatenate(['R0_post'], [str(r) for r in R]))
    # Return to script directory
    os.chdir("..")
    os.chdir("..")
    os.chdir("..")
    os.chdir("scripts")
