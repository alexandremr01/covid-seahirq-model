from testing.utils import Parameters, OutputData
from testing.models import SeahirModel
import numpy as np
from plot_results import plot_fit
import os
import csv

# Teste de compatibilidade para o cenario_generator. O gráfico gerado não deve mudar

I0 = 1
uf = 'SP'
os.chdir('testing')

parameters = Parameters(uf)
model = SeahirModel(uf)
model.run(I0=I0, R0=parameters.r_0, zeta=parameters.zeta, tau=parameters.tau, g=parameters.g, days=parameters.days,
               itv=[10, 9, 8], ifrs=parameters.ifrs, verbose=False)
os.chdir('..')
plot_fit(uf, 150)
exit(0)
