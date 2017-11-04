#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 01:42:29 2017

@author: piotr
"""
import re
import numpy as np
import matplotlib.pyplot as plt


def loadDataFromDataFile(FILE):
    """ """
    dataList = []
    with open(FILE) as infile:
        for line in infile:
            if re.search("#", line):
                continue
            else:
                one = line.split(',')
                ww = re.split(";", one[4])
                dataList.append((float(one[0]), float(one[1]), float(one[2]),
                                 float(one[3]),
                                 np.array(ww[0:-1], dtype=float)))
    return dataList


def getListOfParamsSets(dataList):
    """ """
    paramz = []
    for itm in dataList:
        tupparamz = (itm[0], itm[1])
        if tupparamz not in paramz:
            paramz.append(tupparamz)
    return paramz


def plotOneSet(dataList, paramsTupl, maxxy=(100, 10e5), figNum=1):
    """ """
    ii = 0
    fs = 16
    tkfs = 14
    plt.figure(figNum, figsize=(18, 8))
    plt.subplot(121)
    persist = []
    for itm in dataList:
        tupparamz = (itm[0], itm[1])
        if tupparamz == paramsTupl:
            persist.append(itm[3]/itm[2])
            plt.semilogy(itm[4] + 1.)
            ii += 1
    persist = np.array(persist)
    print("We have", ii, "simulations for", paramsTupl)
    plt.xlim((0, maxxy[0]))
    plt.ylim((0, maxxy[1]))
    plt.xlabel("time [host generations after mutation apperence]", fontsize=fs)
    plt.ylabel("mutant’s relative immunocompetence, $log(y + 1)$ ",
               fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.grid(True)
    plt.subplot(122)
    plt.boxplot(persist)
    plt.ylim((0, 1.0))
    plt.show()
