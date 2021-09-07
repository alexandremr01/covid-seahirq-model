# coding: utf-8
import numpy as np
import csv
import os
import subprocess
from model.scen_gen import scen_gen, TESTING_PARAMETERS_FROM_EPIDEMIOLOGY

class EpidemicPhases:
    def __init__(self, g0, itv0):
        self.g = [g0]
        self.days = [0]
        self.itv_level = [itv0]

    def add_phase(self, g, start_day, itv_level):
        self.g.append(g)
        self.days.append(start_day)
        self.itv_level.append(itv_level)

    def get_length(self):
        return len(self.g)

class SeahirModel:
    """
    Abstrai a chamada do código em C.
    """
    def __init__(self, uf):
        self.uf = uf
        self.cenario_folder = 'cenario' + self.uf        

    def run(self, I0, R0, phases, xI=None, xA=None, testing_parameters=TESTING_PARAMETERS_FROM_EPIDEMIOLOGY, verbose=False, attack=None):
        output_file = 'result_data_' + self.uf + '.csv'
        scen_gen(input_folder=self.cenario_folder, phases=phases, model=3, I_0=I0, R0=R0, fatality=1.0, testing_parameters_mode=testing_parameters, 
                    input_xi=xI, input_xa=xA, verbose=verbose, attack=attack)
        os.chdir("../..")

        if verbose:
            stdout=None
        else:
            stdout=open(os.devnull, 'wb')

        subprocess.call(['bin/csv_to_input', self.cenario_folder], stdout=stdout)
        input_file = 'input/cenarios/' + self.cenario_folder + '/generated-input.txt'
        output_file = '/'.join(['output', output_file])
        subprocess.call(['bin/spatial_covid0d_estrat', input_file, output_file, '3'], stdout=stdout)
        os.chdir("scripts/testing")