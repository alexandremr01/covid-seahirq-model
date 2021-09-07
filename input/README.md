# Input de dados para o executável

O arquivo `template.txt` contém um template para os arquivos de input.

**Atualmente os inputs seguem o formato:**
```
POPULATION_S    EXPOSED_E   INFECTED_I  RECOVERED_R
LAMBDA          BETA        MORT_EQ     MORT_COV        GAMA        A
DAY             LAMBDA      BETA        MORT_EQ         MORT_COV    GAMA       A
DAY             LAMBDA      BETA        MORT_EQ         MORT_COV    GAMA       A
DAY             LAMBDA      BETA        MORT_EQ         MORT_COV    GAMA       A
DAY             LAMBDA      BETA        MORT_EQ         MORT_COV    GAMA       A
...
```

Sendo a primeira linha referênte aos dados iniciais e a segunda linha referente aos parâmetros abaixo:

| Parâmetro     | Significado                           |
|:-------------:|:-------------------------------------:|
| LAMBDA        | Taxa de Natalidade                    |
| BETA          | Taxa de Contaminação                  |
| Mu_equilibrio | Taxa de Mortalidade de Equilibrio     |
| Mu_covid      | Taxa de Mortalidade do COVID-19       |
| GAMA          | Taxa de Recuperação                   |
| a             | Taxa de Conversão E->I                |

Para os valores variando no tempo, valores negativos representam "manter o último valor utilizado"

**Exemplo (template.txt):**
```
# ----------------------------------------------------------
#                   INITIAL VALUES
# ----------------------------------------------------------
# POPULATION_S  EXPOSED_E   INFECTED_I  RECOVERED_R


211049527.0 100.0 1.0 0.0

# ----------------------------------------------------------
#              PARAMETERS (initial value)
# ----------------------------------------------------------
# LAMBDA    BETA    MORT_EQ     MORT_COV    GAMA       A


0.000038 2.950000 0.000018 0.0006 0.016 0.057



# ----------------------------------------------------------
#              PARAMETERS (time changing)
# ----------------------------------------------------------
# DAY   LAMBDA    BETA    MORT_EQ     MORT_COV    GAMA       A


30 -1.00 0.52 -1.00 -1.00 -1.00 -1.00 
31 -1.00 0.5026666667 -1.00 -1.00 -1.00 -1.00 
32 -1.00 0.4691555556 -1.00 -1.00 -1.00 -1.00 
33 -1.00 0.42224 -1.00 -1.00 -1.00 -1.00 
34 -1.00 0.3659413333 -1.00 -1.00 -1.00 -1.00 
35 -1.00 0.3049511111 -1.00 -1.00 -1.00 -1.00 
36 -1.00 0.2439608889 -1.00 -1.00 -1.00 -1.00 
37 -1.00 0.1870366815 -1.00 -1.00 -1.00 -1.00 
38 -1.00 0.1371602331 -1.00 -1.00 -1.00 -1.00 
39 -1.00 0.09601216316 -1.00 -1.00 -1.00 -1.00 
40 -1.00 0.06400810877 -1.00 -1.00 -1.00 -1.00 
41 -1.00 0.04053846889 -1.00 -1.00 -1.00 -1.00 
42 -1.00 0.02432308133 -1.00 -1.00 -1.00 -1.00 
43 -1.00 0.01378307942 -1.00 -1.00 -1.00 -1.00 
44 -1.00 0.007350975692 -1.00 -1.00 -1.00 -1.00 
45 -1.00 0.003675487846 -1.00 -1.00 -1.00 -1.00 
46 -1.00 0.001715227661 -1.00 -1.00 -1.00 -1.00 
47 -1.00 0.00074326532 -1.00 -1.00 -1.00 -1.00 
48 -1.00 0.000297306128 -1.00 -1.00 -1.00 -1.00 
49 -1.00 0.0001090122469 -1.00 -1.00 -1.00 -1.00 
50 -1.00 0.00003633741564 -1.00 -1.00 -1.00 -1.00 
51 -1.00 0.00001090122469 -1.00 -1.00 -1.00 -1.00 
52 -1.00 0.000002906993251 -1.00 -1.00 -1.00 -1.00 
53 -1.00 0.0000006782984253 -1.00 -1.00 -1.00 -1.00 
54 -1.00 0.0000001356596851 -1.00 -1.00 -1.00 -1.00 
55 -1.00 0.00000002260994751 -1.00 -1.00 -1.00 -1.00 
56 -1.00 0.000000003014659668 -1.00 -1.00 -1.00 -1.00 
57 -1.00 0.0000000003014659668 -1.00 -1.00 -1.00 -1.00 
58 -1.00 0 -1.00 -1.00 -1.00 -1.00 
59 -1.00 0 -1.00 -1.00 -1.00 -1.00 
```