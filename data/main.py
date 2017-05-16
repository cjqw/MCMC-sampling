#!env /usr/bin/python3

from math import *

FILE_NAME = 'mh.output'

def parseItem(item):
    x,y = item.split(' ')
    return float(x),float(y)

def parseData(data):
    data = data[2:-2].split('] [')
    return list(map(parseItem,data))

def density(position):
    x,y = position
    f1 = (0.5/pi)*exp(-0.5 * (x * x + y * y))
    x,y = x-1,y-1
    f2 = (0.5/pi)*exp(-0.5 * (x * x + y * y))
    return (f1+f2)/2

def readData():
    with open(FILE_NAME,'r')as fin:
        data = fin.read()
    return parseData(data)

def normalize(item):
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
        index = normalize(item)
        cnt[index] = cnt.get(index,0) + 1
    result = 0
    for key,value in cnt.items():
        p = value/length
        q = density(decode(key)) * 0.01 * 0.01
        result += p*log(p/q)
    return result

data = readData()
for i in range(100):
    length = (i+1) * 10000
    print(length)
    print(KL_Divergence(data[:length]))
