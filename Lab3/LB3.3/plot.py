import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

lines = []
try:
    with open('points.txt') as f:
        lines = f.readlines()
    X = [float(item) for item in lines[0].split()]
    Fi = [float(item) for item in lines[1].split()]
    F1 = [float(item) for item in lines[2].split()]
    F2 = [float(item) for item in lines[3].split()]

    fig, ax = plt.subplots() 
    ax.plot(X, Fi, 'o')
    ax.plot(X, F1,label='1-я степень')
    ax.plot(X, F2,label='2-я степень')
    ax.legend()
    plt.show()
except:
    print('Can\'t open file. Run C++ program to generate.')

