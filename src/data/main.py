#!env /usr/bin/python3

from math import *
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D


FILE_NAME = 'mh.output'
DRAW_SURFACE = True
DRAW_KL = True

def parseItem(item):
    x,y = item.split(' ')
    return float(x),float(y)

def parseData(data):
    data = data[2:-2].split('] [')
    return list(map(parseItem,data))

def density(position):
    x,y = position
    f1 = (0.5/pi)*exp(-0.5 * (x * x + y * y))
    x,y = x-2,y-2
    f2 = (0.5/pi)*exp(-0.5 * (x * x + y * y))
    return (f1+f2)/2

def readData():
    with open(FILE_NAME,'r')as fin:
        data = fin.read()
    return parseData(data)

def encode(item):
    x,y = item
    x = x + 100
    y = y + 100
    x = int(x*100)
    y = int(y*100)
    return x * 1000000 + y

def decode(key):
    x = key//1000000
    y = key - x * 1000000
    x/=100
    y/=100
    return x - 100, y - 100

def KL_Divergence(data):
    length = len(data)
    cnt = {}
    for item in data:
        index = encode(item)
        cnt[index] = cnt.get(index,0) + 1
    result = 0
    for key,value in cnt.items():
        p = value/length
        q = density(decode(key)) * 0.01 * 0.01
        result += p*log(p/q)
    return result

def drawSurface():
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.arange(-3, 6, 0.25)
    Y = np.arange(-3, 6, 0.25)
    X, Y = np.meshgrid(X, Y)
    Z = [[density((x,y)) for x,y in zip(xi,yi)]
         for xi,yi in zip(X,Y)]
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def drawKL(KLValues):
    X,Y = [],[]
    for item in KLValues:
        x,y = item
        X.append(x)
        Y.append(y)
    plt.plot(X,Y,label = "KL-divergence")
    plt.legend()
    plt.show()

if DRAW_SURFACE:
    drawSurface()

data = readData()

KLValues = []
for i in range(0,100):
    length = (i+1) * 10000
    KL_Value = KL_Divergence(data[:length])
    print(length,KL_Value)
    KLValues.append([length,KL_Value])

if DRAW_KL:
    drawKL(KLValues)
# 1000000 0.780989858100924
