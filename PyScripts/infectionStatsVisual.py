#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 01:42:29 2017

@author: piotr
"""
import re
import numpy as np
import matplotlib.pyplot as plt
import infectionStatsBig as isb


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
    allParamsTupl = []
    for itm in dataList:
        tupparamz = (itm[0], itm[1])
        if tupparamz not in allParamsTupl:
            allParamsTupl.append(tupparamz)
    return allParamsTupl


def plotOneSet(dataList, paramsTupl, maxxy=(100, 1e5), figNum=1):
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


def getAvgVals(dataList, allParamsTupl, avgWay='median'):
    """ """
    avgDataList = []
    for prm in allParamsTupl:
        oneGoSimuls = []
        persist = []
        for itm in dataList:
            tupparamz = (itm[0], itm[1])
            if tupparamz == prm:
                persist.append(itm[3]/itm[2])
                oneGoSimuls.append(itm[4])
        meanStats = isb.calcAverageForOneRun(oneGoSimuls, avgWay)
        avgDataList.append((prm[0], prm[1], meanStats, np.array(persist)))
    return sorted(avgDataList, key=lambda elm: (elm[0], elm[1]))


def plotAvgStats(avgDataList, maxxy=(100, 1e5)):
    """ """
    ii = 0
    fs = 16
    tkfs = 14
    linez = ('b-', 'b--', 'r-', 'r--')
    forBoxPlt = []
    boxLbls = []
    plt.figure(1, figsize=(18, 8))
    plt.subplot(121)
    for itm in avgDataList:
        plt.semilogy(itm[2] + 1., linez[ii])
        print("line", linez[ii], "represents", itm[0], ",", itm[1])
        forBoxPlt.append(itm[3])
        boxLbls.append("alpha %1.2f\npath. spp %d" % (itm[0], itm[1]))
        ii += 1
    plt.xlabel("time [host generations after mutation apperence]", fontsize=fs)
    plt.ylabel("mutant’s relative immunocompetence, $log(y + 1)$ ",
               fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.xlim((0, maxxy[0]))
    plt.ylim((0, maxxy[1]))
    plt.grid(True)
    plt.subplot(122)
    boxprops = dict(linestyle='-', linewidth=2.5, color='k')
    medianprops = dict(linestyle='-', linewidth=2.5)
    whiskerprops = dict(linewidth=2.5)
    capprops = dict(linewidth=2.5)
    flierprops = dict(markersize=10)
    plt.boxplot(forBoxPlt, labels=boxLbls, boxprops=boxprops,
                medianprops=medianprops, whiskerprops=whiskerprops,
                capprops=capprops,  flierprops=flierprops)
    plt.ylim((0.5, 1))
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.grid(axis='y')
    plt.show()
