import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Set here
I = [10, 100, 1000, 10000]
show = False

for I0 in I:
    dataframe = pd.read_csv("../../output/testing/xI_table_SP_"+str(I0)+".csv")
    dataframe = dataframe.dropna()

    x = dataframe.xI
    y = dataframe.mortos

    stats = linregress(x, y)

    m = stats.slope
    b = stats.intercept

    print(y)

    plt.scatter(x, y)
    plt.plot(x, m * x + b, color="red")   # I've added a color argument here
    plt.xlabel('xI')
    plt.ylabel('Mortos')
    plt.title('Mortos por xI para I0='+str(I0))
    plt.savefig('../../output/testing/Mortos por xI para I0='+str(I0)+'.png')

    if show:
        plt.show()
    plt.close()
    x = dataframe.xI
    y = dataframe.maxH

    stats = linregress(x, y)

    m = stats.slope
    b = stats.intercept

    print(y)

    plt.scatter(x, y)
    plt.plot(x, m * x + b, color="red")   # I've added a color argument here

    plt.xlabel('xI')
    plt.ylabel('Mortos')
    plt.title('Maximo de hospitalizados por xI para I0='+str(I0))
    plt.savefig('../../output/testing/Max hospitalizados por xI para I0='+str(I0)+'.png')
    if show:
        plt.show()
    plt.close()
