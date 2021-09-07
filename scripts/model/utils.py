import csv
import numpy as np

def read_matrix(filename, age_strata):
    matrix = np.zeros([age_strata, age_strata], dtype=np.float64)
    with open(filename, "r") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        next(spamreader, None)
        j = 0
        for row in spamreader:
            for i in range(age_strata):
                matrix[j, i] = row[i]
            j = j + 1
    return matrix