# coding: utf-8
import numpy as np
import csv
import os
import subprocess
from scipy import optimize
from scipy import stats
import datetime

tables_folder = '../dataset/'
death_table = tables_folder + 'deaths.csv'
hospitalized_table = tables_folder + 'hospitalized.csv'
hospital_deaths_table = tables_folder + 'hospital_deaths.csv'
consolidated_table = tables_folder + 'consolidated.csv'
notifications_folder = '../input/cenarios/'
age_strata = 16


def get_age_range(str):
    """
    Transforma uma string no formato '05-09' em id da faixa etária, inteiro no intervalo [0, age_strata).
    """
    age_range = 0
    try:
        age_range = min(int(str[0:2]) // 5, 15)
    except:
        age_range = -1  # Se str for 'Unknown'
    return age_range


class StateData:
    def __init__(self, uf, num_points):
        """
        Cria um StateData para o estado em uf (ex: SP, AC), lendo do input.
        """
        self.uf = uf
        self.earliest_date = datetime.datetime.strptime('2020-01-01',
                                                        '%Y-%m-%d').date()  # Controla o dia 0 das matrizes.
        self.CD, self.CDday0 = self.read_table(death_table, name='deaths')
        self.H, self.Hday0 = self.read_table(hospitalized_table, name='hospitalized')
        self.C, self.Cday0 = self.read_table(hospital_deaths_table, name='hospitalized deaths')

    def read_table(self, filename, name=''):
        """
        Lê uma tabela CSV e armazena os dados em uma matriz.
        """
        csv_file = open(filename, mode='r')
        stateCSV = csv.reader(csv_file)
        header = next(stateCSV)
        initial_rows = 4
        initial_date = datetime.datetime.strptime(header[initial_rows], '%Y-%m-%d').date()
        m = len(header) - initial_rows - 1  # 1 coluna final que não representa dados
        data = np.zeros((age_strata, m), dtype=np.int32)
        j = 0
        count = 0
        column_uf = 1
        column_age = [idx for idx, s in enumerate(header) if 'age_group' in s][0]  # Coluna que contem "age_group"
        for row in stateCSV:
            if row[column_uf] == self.uf and row[column_age] != "Unknown":
                j = get_age_range(row[column_age])
                data[j, :] += np.array([0 if x == 'NA' else int(x) for x in row[initial_rows:-1]])  # -1
            count += 1
        csv_file.close()
        dates = np.array([datetime.datetime.strptime(x, '%Y-%m-%d').date() for x in header[initial_rows:-1]])
        data = self.fill_empty_dates(data, dates)
        cum_data = np.where(data != 0, np.cumsum(data, axis=1), data)
        # Para pegar a data inicial
        cum_data_total = np.cumsum(data, axis=1)
        cum_data_sum = np.sum(cum_data_total, axis=0)
        date_0 = self.earliest_date + datetime.timedelta(
            days=int(np.where(cum_data_sum > 0)[0][0]))  # Data com 1a ocorrencia
        if name != '':
            print("Initial date from " + name + ": " + str(date_0))
        return cum_data_total, date_0

    def cost(self, outData, D0=0, num_points=150):
        """
        Função de custo para o output dado em outData. Translada D0 dias.
        """
        initial_day = self.Hday0 - datetime.timedelta(days=int(D0 // 1))
        date_delta = (initial_day - self.earliest_date).days
        C = self.C[:, date_delta:(date_delta + num_points)]
        CD = self.CD[:, date_delta:(date_delta + num_points)]
        H = self.H[:, date_delta:(date_delta + num_points)]
        cDays, hDays, cdDays = C.shape[1], CD.shape[1], H.shape[1]
        f_c, f_h, f_cd = outData.get_data_for_days(cDays, hDays, cdDays)
        sumC = self.replace_repeated(C.sum(axis=0).reshape(1, -1))
        sumH = self.replace_repeated(H.sum(axis=0).reshape(1, -1))
        sumCD = self.replace_repeated(CD.sum(axis=0).reshape(1, -1))
        sumf_C = self.replace_repeated(f_c.sum(axis=0).reshape(1, -1))
        sumf_H = self.replace_repeated(f_h.sum(axis=0).reshape(1, -1))
        sumf_CD = self.replace_repeated(f_cd.sum(axis=0).reshape(1, -1))
        f_chi2 = self.normalized_cost_if(f_c, C, H) + self.normalized_cost_if(f_cd, CD, H) \
                 + self.normalized_cost(f_h, H)
        f_chi2 += self.normalized_cost_if(sumf_C, sumC, sumH) + self.normalized_cost_semigauss(sumf_CD, sumCD) \
                  + self.normalized_cost(sumf_H, sumH)
        #print("{:.4g}".format(f_chi2))
        #if f_chi2 == 0:
        #    f_chi2 = 1E23
        return f_chi2

    def cost_deaths(self, outData, D0=0, num_points=150):
        """
        Função de custo para o output dado em outData. Translada D0 dias.
        """
        initial_day = self.Hday0 - datetime.timedelta(days=int(D0 // 1))
        date_delta = (initial_day - self.earliest_date).days
        C = self.C[:, date_delta:(date_delta + num_points)]
        CD = self.CD[:, date_delta:(date_delta + num_points)]
        H = self.H[:, date_delta:(date_delta + num_points)]
        cDays, hDays, cdDays = C.shape[1], CD.shape[1], H.shape[1]
        f_c, f_h, f_cd = outData.get_data_for_days(cDays, hDays, cdDays)
        sumC = self.replace_repeated(C.sum(axis=0).reshape(1, -1))
        sumH = self.replace_repeated(H.sum(axis=0).reshape(1, -1))
        sumCD = self.replace_repeated(CD.sum(axis=0).reshape(1, -1))
        sumf_C = self.replace_repeated(f_c.sum(axis=0).reshape(1, -1))
        #sumf_H = self.replace_repeated(f_h.sum(axis=0).reshape(1, -1))
        sumf_CD = self.replace_repeated(f_cd.sum(axis=0).reshape(1, -1))
        f_chi2 = self.normalized_cost_if(f_cd, CD, H) + self.normalized_cost_if(f_c, C, H)
        f_chi2 += self.normalized_cost_if(sumf_C, sumC, sumH) + self.normalized_cost_semigauss(sumf_CD, sumCD)
        if f_chi2 == 0:
            f_chi2 = 1E23
        #print("{:.4g}".format(f_chi2))
        return f_chi2

    def cost_deaths_ifr(self, outData, D0=0, num_points=150):
        """
        Função de custo para o output dado em outData. Translada D0 dias.
        """
        initial_day = self.Hday0 - datetime.timedelta(days=int(D0 // 1))
        date_delta = (initial_day - self.earliest_date).days
        C = self.C[:, date_delta:(date_delta + num_points)]
        CD = self.CD[:, date_delta:(date_delta + num_points)]
        H = self.H[:, date_delta:(date_delta + num_points)]
        cDays, hDays, cdDays = C.shape[1], CD.shape[1], H.shape[1]
        f_c, f_h, f_cd = outData.get_data_for_days(cDays, hDays, cdDays)
        sumC = self.replace_repeated(C.sum(axis=0).reshape(1, -1))
        sumH = self.replace_repeated(H.sum(axis=0).reshape(1, -1))
        sumCD = self.replace_repeated(CD.sum(axis=0).reshape(1, -1))
        sumf_C = self.replace_repeated(f_c.sum(axis=0).reshape(1, -1))
        # sumf_H = self.replace_repeated(f_h.sum(axis=0).reshape(1, -1))
        sumf_CD = self.replace_repeated(f_cd.sum(axis=0).reshape(1, -1))
        f_chi2 = self.normalized_cost_if(f_cd, CD, H) + self.normalized_cost_if(f_c, C, H)
        if f_chi2 == 0:
            f_chi2 = 1E23
        # print("{:.4g}".format(f_chi2))
        return f_chi2

    def replace_repeated(self, x):
        y = np.zeros((x.shape))
        y[0] = x[0]
        # print(x)
        for i in range(1, len(x)):
            if x[i] != x[i - 1]:
                y[i] = x[i]
        # print(y)
        return y

    def get_d0_from_notif(self):  # Não está em uso
        """
        Retorna data da primeira morte de acordo com notification.csv.
        """
        notifications_file = notifications_folder + 'cenario' + self.uf + '/notification.csv'
        csv_file = open(notifications_file, mode='r')
        notificationCSV = csv.reader(csv_file)
        data = []
        for row in notificationCSV:
            data.append(np.array(row[1:]))
        infected = [-1 if x == '' else int(x) for x in data[0]]
        days_inf = [-1 if x == '' else int(x) for x in data[1]]
        deaths = [-1 if x == '' else int(x) for x in data[2]]
        days_deaths = [-1 if x == '' else int(x) for x in data[3]]
        day_0 = days_deaths[0]
        index_day_0 = -1
        for index, d in enumerate(days_inf):
            if day_0 < d:
                index_day_0 = index - 1
                break
        d0 = datetime.date(day=int(data[4][0]), month=int(data[5][0]), year=int(data[6][0]))
        inf0 = datetime.date(day=int(data[4][0]), month=int(data[5][0]), year=int(data[6][0])) + datetime.timedelta(
            days=days_inf[0] - 1)
        d0 = d0 + datetime.timedelta(days=day_0 - 1)
        print("Initial infection date from notifications: " + str(inf0))
        print("Initial death date from notifications: " + str(d0))
        return d0

    def get_i0(self):  # Não está em uso
        """
        Retorna número de infectados na data especificada em self.get_d0_from_notif com base no notification.csv.
        """
        notifications_file = notifications_folder + 'cenario' + self.uf + '/notification.csv'
        csv_file = open(notifications_file, mode='r')
        notificationCSV = csv.reader(csv_file)
        data = []
        for row in notificationCSV:
            data.append(np.array(row[1:]))
        infected = [-1 if x == '' else int(x) for x in data[0]]
        days_inf = [-1 if x == '' else int(x) for x in data[1]]
        deaths = [-1 if x == '' else int(x) for x in data[2]]
        days_deaths = [-1 if x == '' else int(x) for x in data[3]]
        day_0 = (self.initial_date - datetime.date(day=int(data[4][0]), month=int(data[5][0]),
                                                   year=int(data[6][0]))).days
        index_day_0 = -1
        for index, d in enumerate(days_inf):
            if day_0 < d:
                index_day_0 = index
                break
        infected_day_0 = infected[index_day_0]
        return infected_day_0

    def normalized_cost_semigauss(self, f_y, y):
        array_shape = f_y.shape
        num = np.zeros(array_shape, dtype=np.float64)
        for i in range(0, array_shape[0]):
            for j in range(0, array_shape[1]):
                if f_y[i, j] < y[i, j]:
                    num[i, j] = np.power(f_y[i, j] - y[i, j], 4)
                else:
                    num[i, j] = np.power(f_y[i, j] - y[i, j], 2)
        J = np.divide(num, y, out=np.zeros(num.shape), where=(y != 0))
        return np.nansum(J)

    def normalized_cost(self, f_y, y):
        num = np.power(f_y - y, 2)
        J = np.divide(num, y, out=np.zeros(num.shape), where=(y != 0))
        return np.nansum(J)

    def normalized_cost_if(self, f_y, y, z):
        f_y_sub = f_y[y < z]
        y_sub = y[y < z]
        num = np.power(f_y_sub - y_sub, 2)
        J = np.divide(num, y_sub, out=np.zeros(num.shape), where=(y_sub != 0))
        return np.nansum(J)

    def fill_empty_dates(self, data, dates, n_days=250):
        """
        Preenche dias sem dados com coluna de zeros.
        """
        earliest_date = self.earliest_date
        date_int = np.array([(x - earliest_date).days for x in dates])
        X = np.copy(data[:, 0]).reshape(-1, 1)
        for i in range(1, n_days):
            if i in date_int:
                X = np.append(X, data[:, np.where(date_int == i)[0][0]].reshape(-1, 1), axis=1)
            else:
                X = np.append(X, np.zeros((data.shape[0], 1), dtype=np.int32), axis=1)
        return X

    def get_date_delta(self, D0):
        initial_day = self.Hday0 - datetime.timedelta(days=int(D0 // 1))
        date_delta = (initial_day - self.earliest_date).days
        return date_delta


class OutputFile:
    """
    Lê a saída do modelo especificado em out_file.
    """

    def __init__(self, uf, t_days=400, age_strata=age_strata, compartments=12):
        self.uf = uf
        gamma_h = self.read_gamma_h()
        out_file = '../output/result_data_' + self.uf + '.csv'
        h = np.zeros([t_days, age_strata], dtype=np.float64)
        cd = np.zeros([t_days, age_strata], dtype=np.float64)
        d = np.zeros([t_days, age_strata], dtype=np.float64)
        c = np.zeros([t_days, age_strata], dtype=np.float64)
        y = np.zeros([t_days, age_strata], dtype=np.float64)
        ri = np.zeros([t_days, age_strata], dtype=np.float64)
        c_column = 6
        d_column = 7
        i_column = 2
        h_column = 8
        with open(out_file, "r") as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            next(spamreader, None)
            j = 0
            for row in spamreader:
                for i in range(age_strata):
                    d[j, i] = self.process(row[compartments * (i + 1) + d_column])
                    c[j, i] = self.process(row[compartments * (i + 1) + c_column])
                    y[j, i] = self.process(row[compartments * (i + 1) + i_column])
                    h[j, i] = self.process(row[compartments * (i + 1) + h_column])
                    cd[j, i] = d[j, i] + c[j, i]
                for ii in range(age_strata):
                    ri[j, ii] = self.process(row[compartments * (age_strata + 1) + ii + 1])
                j = j + 1
            # self.C = np.add.accumulate(c, 1)
            # self.CD = np.add.accumulate(cd, 1)
            # self.H = np.add.accumulate(h, 1)
            # print(gamma_h.shape)
            # print(y.shape)
            self.H = np.add.accumulate(np.multiply(gamma_h, y), 0)
            self.CD = np.copy(cd)
            self.C = np.copy(c)

    def read_gamma_h(self):
        param_file = '../input/cenarios/cenario' + self.uf + '/parameters.csv'
        with open(param_file, "r") as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                if len(row) > 0 and row[0] == "GAMA_H":
                    gamma_h = np.array(row[1:], dtype=np.float)
        # print(gamma_h)
        return gamma_h.reshape(1, -1)

    def get_data_for_days(self, cDays, hDays, cdDays):
        return self.C[0:cDays, :].T, self.H[0:hDays, :].T, self.CD[0:cdDays, :].T

    def process(self, x):
        return np.nan if x == "-nan(ind)" else float(x)


def routine_over(x, x_ref, zeta, last_day_out, num_intervention_phases, intervention_levels, stateData):
    # Parametros: R0, tau, zeta, D0, d, g,
    r_0 = x[0]
    tau = x[1]
    ifr_j = x_ref[:]
    n = num_intervention_phases
    g = x[4: n + 4]
    durations = x[n + 4:]
    #print(x)
    days = np.add.accumulate(durations)
    cenario_folder = 'cenario' + stateData.uf
    output_file = 'result_data_over_' + stateData.uf + '.csv'
    I0 = 1
    days = [str(d) for d in days]
    itv = [str(level) for level in intervention_levels]
    # Caso pessimista - cenário selvagem pós flexibilização
    subprocess.call(['python2.7', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', days[0], days[1], days[2],
                     '-m', '3', '-I0', str(I0), '-R0', str(r_0), '-z', str(zeta), '-t', str(tau), '-p', '1',
                     str(g[0]), str(g[1]), str(g[2]), '-itv', '0', itv[0], itv[1], itv[2],
                     '-ex', '0', str(last_day_out), '-ifrs',
                    str(ifr_j[0]), str(ifr_j[1]), str(ifr_j[2]), str(ifr_j[3]), str(ifr_j[4]), str(ifr_j[5]),
                    str(ifr_j[6]), str(ifr_j[7]), str(ifr_j[8]), str(ifr_j[9]), str(ifr_j[10]), str(ifr_j[11]),
                    str(ifr_j[12]), str(ifr_j[13]), str(ifr_j[14]), str(ifr_j[15])], stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(
        ['bin/spatial_covid0d_estrat', 'input/cenarios/' + cenario_folder + '/generated-input.txt',
         '/'.join(['output', output_file]), '3'], stdout=open(os.devnull, 'w'))
    os.chdir("scripts")


def f(x, num_intervention_phases, intervention_levels, stateData, verbose):
    # Parametros: R0, tau, zeta, D0, d, g,
    r_0 = x[0]
    tau = x[1]
    zeta = x[2]
    D0 = x[3]
    n = num_intervention_phases
    g = x[4: n + 4]
    durations = x[n + 4:]
    #print(x)
    days = np.add.accumulate(durations)
    cenario_folder = 'cenario' + stateData.uf
    output_file = 'result_data_' + stateData.uf + '.csv'
    I0 = 1
    days = [str(d) for d in days]
    itv = [str(level) for level in intervention_levels]
    # A proxima linha supõe que num_phases = 3. Precisa ser alterado.
    if verbose:
        subprocess.call(['python2.7', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', days[0], days[1], days[2],
                         '-m', '3', '-I0', str(I0), '-R0', str(r_0), '-z', str(zeta), '-t', str(tau), '-p', '1',
                         str(g[0]),
                         str(g[1]), str(g[2]), '-itv', '0', itv[0], itv[1], itv[2]])  # stdout=open(os.devnull, 'wb'))
    else:
        subprocess.call(['python2.7', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', days[0], days[1], days[2],
                         '-m', '3', '-I0', str(I0), '-R0', str(r_0), '-z', str(zeta), '-t', str(tau), '-p', '1',
                         str(g[0]), str(g[1]), str(g[2]), '-itv', '0', itv[0], itv[1], itv[2]], stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(
        ['bin/spatial_covid0d_estrat', 'input/cenarios/' + cenario_folder + '/generated-input.txt',
         '/'.join(['output', output_file]), '3'], stdout=open(os.devnull, 'w'))
    os.chdir("scripts")
    outData = OutputFile(stateData.uf)
    f_chi2 = stateData.cost(outData, D0)
    #print("{:.4g}".format(f_chi2))
    return f_chi2


def f_ref(x, x_prior, num_intervention_phases, intervention_levels, stateData):
    # Parametros: R0, tau, zeta, D0, d, g,
    r_0 = x_prior[0]
    tau = x_prior[1]
    zeta = x[0]
    D0 = x_prior[3]
    n = num_intervention_phases
    g = x_prior[4: n + 4]
    durations = x_prior[n + 4:]
    ifr_j = x[1:]
    #print(x)
    days = np.add.accumulate(durations)
    cenario_folder = 'cenario' + stateData.uf
    output_file = 'result_data_' + stateData.uf + '.csv'
    I0 = 1
    days = [str(d) for d in days]
    itv = [str(level) for level in intervention_levels]
    subprocess.call(
        ['python', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', days[0], days[1], days[2],
         '-m', '3', '-I0', str(I0), '-R0', str(r_0), '-z', str(zeta), '-t', str(tau), '-p', '1',
         str(g[0]), str(g[1]), str(g[2]), '-itv', '0', itv[0], itv[1], itv[2], '-ifrs',
         str(ifr_j[0]), str(ifr_j[1]), str(ifr_j[2]), str(ifr_j[3]), str(ifr_j[4]), str(ifr_j[5]),
         str(ifr_j[6]), str(ifr_j[7]), str(ifr_j[8]), str(ifr_j[9]), str(ifr_j[10]), str(ifr_j[11]),
         str(ifr_j[12]), str(ifr_j[13]), str(ifr_j[14]), str(ifr_j[15])], stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(
        ['bin/spatial_covid0d_estrat', 'input/cenarios/' + cenario_folder + '/generated-input.txt',
         '/'.join(['output', output_file]), '3'], stdout=open(os.devnull, 'w'))
    os.chdir("scripts")
    outData = OutputFile(stateData.uf)
    f_chi2 = stateData.cost_deaths_ifr(outData, D0)
    # print(f_chi2)
    return f_chi2

def f_ref_ifr(x, x_prior, zeta, num_intervention_phases, intervention_levels, stateData):
    # Parametros: R0, tau, zeta, D0, d, g,
    r_0 = x_prior[0]
    tau = x_prior[1]
    D0 = x_prior[3]
    n = num_intervention_phases
    g = x_prior[4: n + 4]
    durations = x_prior[n + 4:]
    ifr_j = x[:]
    #print(x)
    days = np.add.accumulate(durations)
    cenario_folder = 'cenario' + stateData.uf
    output_file = 'result_data_' + stateData.uf + '.csv'
    I0 = 1
    days = [str(d) for d in days]
    itv = [str(level) for level in intervention_levels]
    subprocess.call(
        ['python', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', days[0], days[1], days[2],
         '-m', '3', '-I0', str(I0), '-R0', str(r_0), '-z', str(zeta), '-t', str(tau), '-p', '1',
         str(g[0]), str(g[1]), str(g[2]), '-itv', '0', itv[0], itv[1], itv[2], '-ifrs',
         str(ifr_j[0]), str(ifr_j[1]), str(ifr_j[2]), str(ifr_j[3]), str(ifr_j[4]), str(ifr_j[5]),
         str(ifr_j[6]), str(ifr_j[7]), str(ifr_j[8]), str(ifr_j[9]), str(ifr_j[10]), str(ifr_j[11]),
         str(ifr_j[12]), str(ifr_j[13]), str(ifr_j[14]), str(ifr_j[15])], stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(
        ['bin/spatial_covid0d_estrat', 'input/cenarios/' + cenario_folder + '/generated-input.txt',
         '/'.join(['output', output_file]), '3'], stdout=open(os.devnull, 'w'))
    os.chdir("scripts")
    outData = OutputFile(stateData.uf)
    f_chi2 = stateData.cost_deaths(outData, D0)
    # print(f_chi2)
    return f_chi2


class Optimizer:
    def __init__(self, verbose=True):
        self.verbose = verbose

    def run(self, uf, num_intervention_phases=3, intervention_levels=[10, 9, 8], num_points=150, n_iteration=20, pop_size=30):
        r_bound = [2.0, 4.0]
        g_bound = [0.01, 1.0]
        tau_bond = self.get_tau_bounds(uf)
        print("Tau bonds:")
        print(tau_bond)
        day_bound = [3.0, 60.0]
        zeta_bound = [0.01, 0.8]
        d0_bound = [3, 20.0]
        bounds = [r_bound, tau_bond, zeta_bound, d0_bound]
        for i in range(num_intervention_phases):
            bounds.append(g_bound)
        for i in range(num_intervention_phases):
            bounds.append(day_bound)
        print("Otimização para " + uf)
        stateData = StateData(uf, num_points)
        ret = optimize.differential_evolution(f, bounds, args=(
            num_intervention_phases, intervention_levels, stateData, self.verbose), disp=True,
                                              maxiter=n_iteration, popsize=pop_size, seed=1234)
        print(ret.x)
        # Garante que o result_data_UF.csv foi gerado com os parâmetros mais atualizados
        final_cost = f(ret.x, num_intervention_phases, intervention_levels, stateData, self.verbose)
        print("Valor final: %f" % final_cost)
        self.save_parameters('cenario' + uf, ret.x, stateData)
        return ret

    def run_tune(self, uf, ret, num_intervention_phases=3, intervention_levels=[10, 9, 8], num_points=150,
                 n_iteration_ref=20, pop_size=30):
        zeta_bound = [0.01, 0.8]
        bounds = [zeta_bound]
        ifr_bound = [[0.00001, 0.0004], [0.00001, 0.0004], [0.00001, 0.0005], [0.00001, 0.0002], [0.00003, 0.002],
                     [0.00005, 0.002],
                     [0.0001, 0.008], [0.0001, 0.008], [0.0002, 0.01], [0.0004, 0.01], [0.001, 0.05], [0.001, 0.05],
                     [0.005, 0.05],
                     [0.01, 0.1], [0.02, 0.2], [0.05, 0.2]]
        for i in range(0, age_strata):
            bounds.append(ifr_bound[i])
        print("Refinamento da otimização para " + uf)
        stateData = StateData(uf, num_points)
        ret_ref = optimize.differential_evolution(f_ref, bounds, args=(
            ret.x, num_intervention_phases, intervention_levels, stateData), disp=True,
                                                  maxiter=n_iteration_ref, popsize=pop_size, seed=1234)
        print(ret_ref.x)
        # Garante que o result_data_UF.csv foi gerado com os parâmetros mais atualizados
        final_cost = f_ref(ret_ref.x, ret.x, num_intervention_phases, intervention_levels, stateData)
        print("Valor final: %f" % final_cost)
        return ret_ref

    def run_tune_ifr(self, uf, ret, zeta, num_intervention_phases=3, intervention_levels=[10, 9, 8], num_points=150,
                 n_iteration_ref=20, pop_size=30):
        bounds = [[0.00001, 0.0004], [0.00001, 0.0004], [0.00001, 0.0005], [0.00001, 0.0002], [0.00003, 0.002],
                     [0.00005, 0.002],
                     [0.0001, 0.008], [0.0001, 0.008], [0.0002, 0.01], [0.0004, 0.01], [0.001, 0.05], [0.001, 0.05],
                     [0.005, 0.05],
                     [0.01, 0.1], [0.02, 0.2], [0.05, 0.2]]
        print("Refinamento dos IFRs para " + uf)
        stateData = StateData(uf, num_points)
        ret_ref_ifr = optimize.differential_evolution(f_ref_ifr, bounds, args=(
            ret.x, zeta, num_intervention_phases, intervention_levels, stateData), disp=True,
                                                  maxiter=n_iteration_ref, popsize=pop_size, seed=1234)
        print(ret_ref_ifr.x)
        # Garante que o result_data_UF.csv foi gerado com os parâmetros mais atualizados
        final_cost = f_ref_ifr(ret_ref_ifr.x, ret.x, zeta, num_intervention_phases, intervention_levels, stateData)
        print("Valor final: %f" % final_cost)
        # Chama rotina para armazenar resultados em um cenário pessimista
        last_day_out = stateData.get_date_delta(ret.x[3]) + num_points
        routine_over(ret.x, ret_ref_ifr.x, zeta, last_day_out, num_intervention_phases, intervention_levels, stateData)
        self.save_parameters_ref('cenario' + uf, ret.x, ret_ref_ifr.x, zeta, stateData)
        return ret_ref_ifr

    def get_tau_bounds(self, uf):
        epidemiologic_file = '../input/cenarios/cenario' + uf + '/epidemiology_data.csv'
        with open(epidemiologic_file, "r") as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            epid_variables = 12
            age_strata = 16
            next(spamreader, None)
            j = 0
            data_epid = np.zeros([epid_variables, age_strata], dtype=np.float64)
            for row in spamreader:
                for i in range(age_strata):
                    data_epid[j, i] = row[i + 1]
                j = j + 1
        TC = np.max(data_epid[0, :])
        upper_bound = 1.0 / TC
        return [1.1, 0.9 * upper_bound]

    def save_parameters(self, cenario_folder, x, stateData):
        optimized_parameters_csv = os.path.join('..', 'input', 'cenarios', cenario_folder, 'optimized_parameters.csv')
        with open(optimized_parameters_csv, 'w') as csv_file:
            spamwriter = csv.writer(csv_file)
            spamwriter.writerow(np.concatenate((['r_0'], [str(x[0])])))
            spamwriter.writerow(np.concatenate((['tau'], [str(x[1])])))
            spamwriter.writerow(np.concatenate((['zeta'], [str(x[2])])))
            spamwriter.writerow(np.concatenate((['D0'], [str(x[3])])))
            spamwriter.writerow(np.concatenate((['g1'], [str(x[4])])))
            spamwriter.writerow(np.concatenate((['g2'], [str(x[5])])))
            spamwriter.writerow(np.concatenate((['g3'], [str(x[6])])))
            spamwriter.writerow(np.concatenate((['d1'], [str(x[7])])))
            spamwriter.writerow(np.concatenate((['d2'], [str(x[8])])))
            spamwriter.writerow(np.concatenate((['d3'], [str(x[9])])))
            # spamwriter.writerow(np.concatenate((['initial_optimization_date'], [str(stateData.initial_date)])))

    def save_parameters_ref(self, cenario_folder, x, x_ref, zeta, stateData):
        optimized_parameters_csv = os.path.join('..', 'input', 'cenarios', cenario_folder,
                                                'optimized_parameters_tuned.csv')
        with open(optimized_parameters_csv, 'w') as csv_file:
            spamwriter = csv.writer(csv_file)
            spamwriter.writerow(np.concatenate((['r_0'], [str(x[0])])))
            spamwriter.writerow(np.concatenate((['tau'], [str(x[1])])))
            spamwriter.writerow(np.concatenate((['zeta'], [zeta])))
            spamwriter.writerow(np.concatenate((['D0'], [str(x[3])])))
            spamwriter.writerow(np.concatenate((['g1'], [str(x[4])])))
            spamwriter.writerow(np.concatenate((['g2'], [str(x[5])])))
            spamwriter.writerow(np.concatenate((['g3'], [str(x[6])])))
            spamwriter.writerow(np.concatenate((['d1'], [str(x[7])])))
            spamwriter.writerow(np.concatenate((['d2'], [str(x[8])])))
            spamwriter.writerow(np.concatenate((['d3'], [str(x[9])])))
            spamwriter.writerow(np.concatenate((['ifr_1'], [str(x_ref[0])])))
            spamwriter.writerow(np.concatenate((['ifr_2'], [str(x_ref[1])])))
            spamwriter.writerow(np.concatenate((['ifr_3'], [str(x_ref[2])])))
            spamwriter.writerow(np.concatenate((['ifr_4'], [str(x_ref[3])])))
            spamwriter.writerow(np.concatenate((['ifr_5'], [str(x_ref[4])])))
            spamwriter.writerow(np.concatenate((['ifr_6'], [str(x_ref[5])])))
            spamwriter.writerow(np.concatenate((['ifr_7'], [str(x_ref[6])])))
            spamwriter.writerow(np.concatenate((['ifr_8'], [str(x_ref[7])])))
            spamwriter.writerow(np.concatenate((['ifr_9'], [str(x_ref[8])])))
            spamwriter.writerow(np.concatenate((['ifr_10'], [str(x_ref[9])])))
            spamwriter.writerow(np.concatenate((['ifr_11'], [str(x_ref[10])])))
            spamwriter.writerow(np.concatenate((['ifr_12'], [str(x_ref[11])])))
            spamwriter.writerow(np.concatenate((['ifr_13'], [str(x_ref[12])])))
            spamwriter.writerow(np.concatenate((['ifr_14'], [str(x_ref[13])])))
            spamwriter.writerow(np.concatenate((['ifr_15'], [str(x_ref[14])])))
            spamwriter.writerow(np.concatenate((['ifr_16'], [str(x_ref[15])])))
