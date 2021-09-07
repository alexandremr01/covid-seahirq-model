# coding: utf-8
import numpy as np
import csv
import os
import datetime

# Global parameters
tables_folder = '../dataset/'
death_table = tables_folder + 'deaths.csv'
hospitalized_table = tables_folder + 'hospitalized.csv'
hospital_deaths_table = tables_folder + 'hospital_deaths.csv'
consolidated_table = tables_folder + 'consolidated.csv'
notifications_folder = '../input/cenarios/'
age_strata = 16
t_days = 400
compartments = 12


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
    def __init__(self, uf):
        """
        Cria um StateData para o estado em uf (ex: SP, AC), lendo do input.
        """
        self.uf = uf
        self.earliest_date = datetime.datetime.strptime('2020-01-01',
                                                        '%Y-%m-%d').date()  # Controla o dia 0 das matrizes.
        self.CD, self.CDday0, self.CDLastDay = self.read_table(death_table, name='deaths')
        self.H, self.Hday0, self.HLastDay = self.read_table(hospitalized_table, name='hospitalized')
        self.C, self.Cday0, self.CLastDay = self.read_table(hospital_deaths_table, name='hospitalized deaths')
        # print(np.sum(self.CD[:,], axis=0))
        self.d_last = self.get_last_day_deaths()

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
                data[j, :] += np.array([0 if x == 'NA' else int(x) for x in row[initial_rows:-1]])
            count += 1
        csv_file.close()
        dates = np.array([datetime.datetime.strptime(x, '%Y-%m-%d').date() for x in header[initial_rows:-1]])
        last_date = dates[-1]
        data = self.fill_empty_dates(data, dates)
        cum_data = np.cumsum(data, axis=1)
        # Para pegar a data inicial
        cum_data_sum = np.sum(cum_data, axis=0)
        # index_date_0 = np.where(cum_data_sum>0)[0][0] + initial_rows  # Data com 1a ocorrencia
        # date_0 = datetime.datetime.strptime(header[index_date_0], '%Y-%m-%d').date()
        date_0 = self.earliest_date + datetime.timedelta(days=int(np.where(cum_data_sum > 0)[0][0]))
        if name != '':
            print("Initial date from " + name + ": " + str(date_0))
        return cum_data, date_0, last_date

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

    def trim(self, parameters):
        self.initial_day = self.Hday0 - datetime.timedelta(days=int(parameters.D0 // 1))
        date_delta = (self.initial_day - self.earliest_date).days
        self.C = np.copy(self.C[:, date_delta:])
        self.CD = np.copy(self.CD[:, date_delta:])
        self.H = np.copy(self.H[:, date_delta:])

    def get_last_day_deaths(self):
        date_delta = (self.CDLastDay - self.earliest_date).days
        return np.sum(self.CD[:, date_delta])


class Parameters:
    """
    Lê os parâmetros salvos.
    """

    def __init__(self, uf, num_intervention_phases=3):
        self.uf = uf
        cenario = '../input/cenarios/cenario' + uf + '/optimized_parameters.csv'
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
        durations = [float(x) for x in parameters[n + 4:]]
        self.days = np.add.accumulate(durations)
        self.gamma_h = self.read_gamma_h()

    def read_gamma_h(self):
        param_file = '../input/cenarios/cenario' + self.uf + '/parameters.csv'
        with open(param_file, "r") as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                if len(row) > 0 and row[0] == "GAMA_H":
                    gamma_h = np.array(row[1:], dtype=np.float)
        # print(gamma_h)
        return gamma_h.reshape(1, -1)


class ResultsData:
    """
    Lê o output do modelo.
    """

    def __init__(self, uf, param, flag=1):
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
        self.t_array = np.array(t_days, dtype=np.float64)

        self.t_array = np.linspace(0, t_days - 1, t_days)

        if flag == 1:
            output_file = '../output/result_data_' + uf + '.csv'
        else:
            output_file = '../output/result_data_over_' + uf + '.csv'

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
        print(param.gamma_h.shape)
        print(self.y.shape)
        self.hac = np.add.accumulate(np.multiply(param.gamma_h, self.y), axis=0)
        self.HAc = np.sum(self.hac, axis=1)
        self.L = np.sum(self.l, axis=1)
        self.RI = np.sum(self.ri, axis=1)
        self.CD = self.C + self.D
        # print(self.CD)
        self.Ac = self.C + self.RI + self.I + self.D


def save_extrapolation(cenario_folder, x, sum_x):
    optimized_parameters_csv = os.path.join('..', 'input', 'cenarios', cenario_folder, 'extrapolated_parameters.csv')
    with open(optimized_parameters_csv, 'w') as csv_file:
        spamwriter = csv.writer(csv_file)
        spamwriter.writerow(np.concatenate((['Age range'], ['Model last deaths'], ['Data last death'],
                                            ['Relative error'], ['Death last day/year'],
                                            ['Death last day/year - Pessimistic'])))
        spamwriter.writerow(np.concatenate(
            (['0 to 4'], [str(x[0, 0])], [str(x[0, 1])], [str(x[0, 2])], [str(x[0, 3])], [str(x[0, 4])])))
        spamwriter.writerow(np.concatenate(
            (['5 to 9'], [str(x[1, 0])], [str(x[1, 1])], [str(x[1, 2])], [str(x[1, 3])], [str(x[1, 4])])))
        spamwriter.writerow(np.concatenate(
            (['10 to 14'], [str(x[2, 0])], [str(x[2, 1])], [str(x[2, 2])], [str(x[2, 3])], [str(x[2, 4])])))
        spamwriter.writerow(np.concatenate(
            (['15 to 19'], [str(x[3, 0])], [str(x[3, 1])], [str(x[3, 2])], [str(x[3, 3])], [str(x[3, 4])])))
        spamwriter.writerow(np.concatenate(
            (['20 to 24'], [str(x[4, 0])], [str(x[4, 1])], [str(x[4, 2])], [str(x[4, 3])], [str(x[4, 4])])))
        spamwriter.writerow(np.concatenate(
            (['25 to 29'], [str(x[5, 0])], [str(x[5, 1])], [str(x[5, 2])], [str(x[5, 3])], [str(x[5, 4])])))
        spamwriter.writerow(np.concatenate(
            (['30 to 34'], [str(x[6, 0])], [str(x[6, 1])], [str(x[6, 2])], [str(x[6, 3])], [str(x[6, 4])])))
        spamwriter.writerow(np.concatenate(
            (['35 to 39'], [str(x[7, 0])], [str(x[7, 1])], [str(x[7, 2])], [str(x[7, 3])], [str(x[7, 4])])))
        spamwriter.writerow(np.concatenate(
            (['40 to 44'], [str(x[8, 0])], [str(x[8, 1])], [str(x[8, 2])], [str(x[8, 3])], [str(x[8, 4])])))
        spamwriter.writerow(np.concatenate(
            (['45 to 49'], [str(x[9, 0])], [str(x[9, 1])], [str(x[9, 2])], [str(x[9, 3])], [str(x[9, 4])])))
        spamwriter.writerow(np.concatenate(
            (['50 to 54'], [str(x[10, 0])], [str(x[10, 1])], [str(x[10, 2])], [str(x[10, 3])], [str(x[10, 4])])))
        spamwriter.writerow(np.concatenate(
            (['55 to 59'], [str(x[11, 0])], [str(x[11, 1])], [str(x[11, 2])], [str(x[11, 3])], [str(x[11, 4])])))
        spamwriter.writerow(np.concatenate(
            (['60 to 64'], [str(x[12, 0])], [str(x[12, 1])], [str(x[12, 2])], [str(x[12, 3])], [str(x[12, 4])])))
        spamwriter.writerow(np.concatenate(
            (['65 to 69'], [str(x[13, 0])], [str(x[13, 1])], [str(x[13, 2])], [str(x[13, 3])], [str(x[13, 4])])))
        spamwriter.writerow(np.concatenate(
            (['70 to 74'], [str(x[14, 0])], [str(x[14, 1])], [str(x[14, 2])], [str(x[14, 3])], [str(x[14, 4])])))
        spamwriter.writerow(np.concatenate(
            (['75 or more'], [str(x[15, 0])], [str(x[15, 1])], [str(x[15, 2])], [str(x[15, 3])], [str(x[15, 4])])))
        spamwriter.writerow(np.concatenate(
            (['Soma'], [str(sum_x[0])], [str(sum_x[1])], [str(sum_x[2])], [str(sum_x[3])], [str(sum_x[4])])))


def extrapolate_var(uf, num_intervention_phases):
    cenario_folder = 'cenario' + uf

    parameters = Parameters(uf, num_intervention_phases)
    resData = ResultsData(uf, parameters, 1)
    resData_pess = ResultsData(uf, parameters, 0)
    stateData = StateData(uf)

    stateData.trim(parameters)
    print("Morte no último dia de dados: %d" % stateData.d_last)
    d_real = stateData.d_last
    last_day_out = (stateData.CDLastDay - stateData.initial_day).days
    last_day_date = datetime.datetime.strptime('2020-12-31', '%Y-%m-%d').date()
    last_day_year = (last_day_date - stateData.initial_day).days
    d_prev = np.around(resData.CD[last_day_out])
    d_last_year = np.around(resData.CD[last_day_year])
    d_last_year_pess = np.around(resData_pess.CD[last_day_year])
    percentual = (d_real - d_prev) / d_real

    percentual_array = np.divide((stateData.CD[:, -1] - resData.cd.T[:, last_day_out]), stateData.CD[:, -1])
    x = np.vstack((np.around(resData.cd.T[:, last_day_out]), stateData.CD[:, -1], percentual_array,
                  np.around(resData.cd.T[:, last_day_year]), np.around(resData_pess.cd.T[:, last_day_year])))
    sum_x = np.concatenate(([d_prev], [d_real], [percentual], [d_last_year], [d_last_year_pess]), axis=0)
    save_extrapolation(cenario_folder, x.T, sum_x)

    print("Morte prevista pelo modelo no último dia de dados: %d" % np.around(resData.CD[last_day_out]))
    print("Erro percentual do modelo no último dia de dados: %f" % percentual)
    print("Morte prevista pelo modelo no último dia do ano: %d" % np.around(resData.CD[last_day_year]))
    print("Morte prevista pelo modelo no último dia do ano (caso pessimiste): %d" % np.around(resData_pess.CD[last_day_year]))
