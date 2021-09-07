# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import datetime

tables_folder = '../dataset/'
death_table = tables_folder + 'deaths.csv'
hospitalized_table = tables_folder + 'hospitalized.csv'
hospital_deaths_table = tables_folder + 'hospital_deaths.csv'
consolidated_table = tables_folder + 'consolidated.csv'
notifications_folder = '../input/cenarios/'
age_strata = 16
t_days = 400
compartments = 12
num_intervention_phases = 3


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


def trim_X(X, day0, date):
    delta = (date - day0).days
    Y = np.copy(X[:, delta:])
    return Y


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

    def trim(self, parameters, num_points):
        self.initial_day = self.Hday0 - datetime.timedelta(days=int(parameters.D0 // 1))
        date_delta = (self.initial_day - self.earliest_date).days
        self.C = np.copy(self.C[:, date_delta:(date_delta+num_points)])
        self.CD = np.copy(self.CD[:, date_delta:(date_delta+num_points)])
        self.H = np.copy(self.H[:, date_delta:(date_delta+num_points)])

    def get_last_day_deaths(self):
        date_delta = (self.CDLastDay - self.earliest_date).days
        return np.sum(self.CD[:, date_delta])


class NotificationData:
    def __init__(self, uf):
        notifications_file = '../input/cenarios/cenario' + uf + '/notification.csv'
        with open(notifications_file, "r") as csv_file:
            notificationCSV = csv.reader(csv_file)
            data = []
            for row in notificationCSV:
                data.append(np.array(row[1:]))
            infected = np.array([-1 if x == '' else int(x) for x in data[0]])
            days_inf = np.array([-1 if x == '' else int(x) for x in data[1]])
            deaths = np.array([-1 if x == '' else int(x) for x in data[2]])
            days_deaths = np.array([-1 if x == '' else int(x) for x in data[3]])
            self.leitos = int(data[7][0])
            self.date0 = datetime.date(day=int(data[4][0]), month=int(data[5][0]), year=int(data[6][0]))
        num_days = days_inf[-1]
        YData = np.zeros((num_days), dtype=np.int)
        days_inf = days_inf - 1
        YData[days_inf] = infected
        self.day_0 = days_deaths[0]
        self.YData = fill_zeros(YData)
        # self.YData = YData[self.day_0-1:]

    def trim(self, parameters, stateData):
        self.initial_day = stateData.Hday0 - datetime.timedelta(days=parameters.D0)
        date_delta = (self.initial_day - self.date0).days
        # print(date_delta)
        if date_delta > 0:
            self.YData = np.copy(self.YData[:, date_delta:])
        else:
            new_ydata = np.append(np.zeros((1, -date_delta), dtype=np.int32), self.YData)
            self.YData = np.copy(new_ydata)


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


class OutputData:
    """
    Lê o output do modelo.
    """

    def __init__(self, uf, param):
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

        output_file = '../output/result_data_' + uf + '.csv'

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


def fill_zeros(X):
    previous = X[0]
    Y = X
    for i, xi in enumerate(X):
        if xi == 0:
            Y[i] = previous
        else:
            previous = xi
    return Y


def save_fig(suffix_fig, uf):
    cenario_folder = 'cenario' + uf
    path_to_figure = '../output/figures'
    dirs = os.path.join(path_to_figure, cenario_folder)
    try:
        os.makedirs(dirs)
    except OSError:
        pass
    plt.savefig(os.path.join(dirs, cenario_folder[:] + suffix_fig), format='svg',
                dpi=300, bbox_inches='tight')
    return 0


def plot_fit(uf, num_points):
    output_file = 'result_data_' + uf + '.csv'

    stateData = StateData(uf)
    notifData = NotificationData(uf)
    parameters = Parameters(uf, num_intervention_phases)
    outData = OutputData(uf, parameters)

    notifData.trim(parameters, stateData)
    stateData.trim(parameters, num_points)
    # print("Morte no último dia de dados: %d"%stateData.d_last)
    d_real = stateData.d_last
    last_day_out = (stateData.CDLastDay - stateData.initial_day).days
    d_prev = outData.CD[last_day_out]
    percentual = (d_real - d_prev) / d_real

    # print("Morte prevista pelo modelo no último dia de dados: %d"%outData.CD[last_day_out])
    # print("Erro percentual do modelo no último dia de dados: %f" % percentual)

    CDData = np.sum(stateData.CD, axis=0)
    CData = np.sum(stateData.C, axis=0)
    HData = np.sum(stateData.H, axis=0)

    tData = np.linspace(0, notifData.YData.size - 1, notifData.YData.size)
    cData = np.linspace(0, CData.size - 1, CData.size)
    hData = np.linspace(0, HData.size - 1, HData.size)

    # Criam-se grupos de faixa etária: 1: 0 - 20 / 2: 20 - 55 / 3: 55 ou mais

    I1 = np.sum(outData.y[:, 0:3], axis=1)
    I2 = np.sum(outData.y[:, 4:9], axis=1)
    I3 = np.sum(outData.y[:, 10:16], axis=1)

    R1 = np.sum(outData.r[:, 0:3], axis=1)
    R2 = np.sum(outData.r[:, 4:9], axis=1)
    R3 = np.sum(outData.r[:, 10:16], axis=1)

    H1 = np.sum(outData.h[:, 0:3], axis=1)
    H2 = np.sum(outData.h[:, 4:9], axis=1)
    H3 = np.sum(outData.h[:, 10:16], axis=1)

    L1 = np.sum(outData.l[:, 0:3], axis=1)
    L2 = np.sum(outData.l[:, 4:9], axis=1)
    L3 = np.sum(outData.l[:, 10:16], axis=1)

    C1 = np.sum(outData.c[:, 0:3], axis=1)
    C2 = np.sum(outData.c[:, 4:9], axis=1)
    C3 = np.sum(outData.c[:, 10:16], axis=1)

    # Manter compatibilidade com versão antiga
    day_init = 0
    day_next_1 = parameters.days[0]
    day_next_2 = parameters.days[1]
    day_next_3 = parameters.days[2]

    # Plota os gráficos
    plt.figure(1)
    plt.subplot(121)
    t_array = outData.t_array
    # plt.plot(t_array, outData.S, '-', label='S')
    plt.plot(t_array, outData.I, '-', label='I')
    # plt.plot(t_array, outData.R, '-', label='R')
    plt.plot(t_array, outData.CD, '-r', label='CD')
    plt.plot(t_array, outData.C, '-g', label='C')
    # plt.plot(t_array, outData.C, '-', label='C')
    plt.plot(t_array, outData.H, '-', label='H')
    plt.plot(t_array, outData.HAc, '-', label='HAc')
    # plt.plot(t_array, outData.N, '-', label='N')
    # plt.plot(t_array, outData.E, '-', label='E')
    # plt.plot(t_array, outData.A, '-', label='A')
    # plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='blue')
    plt.axvspan(day_next_2, day_next_3, alpha=0.1, color='red')
    plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='black')
    plt.xlim([0, 0.7 * t_days])
    plt.xlabel(u'dias')
    plt.ylabel(u'indivíduos')
    plt.legend(loc='center left', shadow=True, fontsize='small')

    plt.subplot(122)
    # plt.semilogy(t_array, outData.S)
    # plt.semilogy(t_array, outData.S)
    plt.semilogy(t_array, outData.I)
    # plt.semilogy(t_array, outData.R)
    plt.semilogy(t_array, outData.CD, '-r')
    plt.semilogy(t_array, outData.C, '-g')
    # plt.semilogy(t_array, outData.C)
    plt.semilogy(t_array, outData.H)
    plt.semilogy(t_array, outData.HAc, '-b')
    plt.semilogy(tData, notifData.YData, 'ok')
    plt.semilogy(cData, CDData, 'or', label=u'Óbitos totais')
    plt.semilogy(hData, HData, 'ob', label=u'Hospitalizados')
    plt.semilogy(cData, CData, 'og', label=u'Óbitos hosp')

    plt.xlim([0, 0.7 * t_days])
    plt.ylim([1, 1.1 * outData.Ac.max()])
    plt.xlabel('tempo (dias)')
    # plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='blue')
    plt.axvspan(day_next_2, day_next_3, alpha=0.1, color='red')
    plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='black')
    max_H = np.max(outData.H)
    max_L = np.max(outData.L)
    max_I = np.max(outData.I)
    max_C = np.max(outData.CD)
    t_max_I = t_array[np.where(outData.I == max_I)]
    t_max_L = t_array[np.where(outData.L == max_L)]
    textstr = '\n'.join([r'$Max(H)=%.2e$' % (max_H,), r'$Max(L)=%.2e$' % (max_L,), r'$Max(I)=%.2e$' % (max_I,),
                         r'$t(max(I))=%.f$ dias' % (t_max_I,),
                         r'$t(max(L))=%.f$ dias' % (t_max_L,),
                         r'Obitos estimados $=%.2e$' % (max_C,),
                         'dia zero: 26/02'])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    plt.text(0.5, 0.2, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)

    plt.suptitle(u'Curva total das populações em compartimentos, ' + uf + ' - modelo SEAHIR')

    #plt.show()

    _, ax = plt.subplots(figsize=(6, 6))
    # plt.semilogy(t_array, outData.S)
    plt.semilogy(t_array, outData.I, label=u'Infectados atuais - modelo')
    # plt.semilogy(t_array, outData.R)
    plt.semilogy(t_array, outData.CD, '-r', label=u'Óbitos totais - modelo')
    plt.semilogy(t_array, outData.C, '-g', label=u'Óbitos hosp - modelo')
    # plt.semilogy(t_array, outData.C)
    plt.semilogy(t_array, outData.H)
    plt.semilogy(t_array, outData.HAc, '-b')
    plt.semilogy(tData, notifData.YData, 'ok')
    plt.semilogy(cData, CDData, 'or', label=u'Óbitos totais')
    plt.semilogy(hData, HData, 'ob', label=u'Hospitalizados')
    plt.semilogy(cData, CData, 'og', label=u'Óbitos hosp')

    plt.xlim([0, 0.7 * t_days])
    plt.ylim([1, 1.1 * outData.Ac.max()])
    plt.xlabel('tempo (dias)')
    # plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='blue')
    plt.axvspan(day_next_2, day_next_3, alpha=0.1, color='red')
    plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='black')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    plt.legend(loc='lower right', shadow=True, fontsize='small')

    plt.suptitle(u'Curva total das populações em compartimentos, ' + uf + ' - modelo SEAHIR')
    save_fig('fit', uf)
    plt.show()