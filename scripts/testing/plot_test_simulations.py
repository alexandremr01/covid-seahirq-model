# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import datetime
from utils import OutputData

t_days=400

def save_fig(suffix_fig, uf):
    cenario_folder = 'cenario' + uf
    path_to_figure = '../output/figures'
    dirs = os.path.join(path_to_figure, cenario_folder)
    try:
        os.makedirs(dirs)
    except OSError:
        pass
    plt.savefig(os.path.join(dirs, cenario_folder[:] + suffix_fig), format='png',
                dpi=300, bbox_inches='tight')
    return 0


def plot_fit(uf, parameters, day_to_start_testing=0, xI=0, text=None):
    output_file = 'result_data_' + uf + '.csv'
    outData = OutputData(uf, parameters)

    # print("Morte prevista pelo modelo no último dia de dados: %d")

    tData = np.linspace(0, t_days-1, t_days)
    cData = np.linspace(0, t_days - 1, t_days)
    hData = np.linspace(0, t_days - 1, t_days)

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
    plt.xlim([0, t_days])
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
    plt.xlim([0, t_days])
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

    # plt.show()

    _, ax = plt.subplots(figsize=(6, 6))
    # plt.semilogy(t_array, outData.S)
    plt.semilogy(t_array, outData.I, label=u'Infectados atuais')
    # plt.semilogy(t_array, outData.R)
    plt.semilogy(t_array, outData.CD, '-r', label=u'Óbitos totais')
    # plt.semilogy(t_array, outData.C, '-g', label=u'Óbitos hosp')
    plt.semilogy(t_array, outData.Qi, label=u'Quarentenados infectados')
    # plt.semilogy(t_array, outData.Qa, label=u'Quarentenados assintomáticos')
    plt.semilogy(t_array, outData.Qi + outData.I, label=u'Infectados - quarentenados ou não')
    # plt.semilogy(t_array, outData.C)
    # plt.semilogy(t_array, outData.H, label=u'Hospitalizados')
    plt.xlim([0, t_days])
    plt.ylim([1, 1.1 * outData.Ac.max()])
    plt.xlabel('tempo (dias)')
    # plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='blue')
    plt.axvspan(day_next_2, day_next_3, alpha=0.1, color='red')
    plt.axvspan(day_next_1, day_next_2, alpha=0.1, color='black')

    plt.axvline(x=day_to_start_testing, ymin=0, ymax= 1.1 * outData.I.max(), color='black')

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    plt.legend(loc='lower right', shadow=True, fontsize='small')
    if text is not None:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(0.5, 0.8, text, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center',
                 bbox=props)

    title = u'Curva total das populações em compartimentos, ' + uf + ' - modelo SEAHIRQ'
    if xI is not None:
        title += '\nxI='+str(xI)
    plt.suptitle(title)
    save_fig('.png', uf)
    plt.show()