#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Your doc string here, please...


Created on Thu Apr  6 16:44:53 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


datype = np.dtype([('VAR', 'f8'), ('Mode', 'S12'), ('meanAllel', 'f8'),
                   ('stdAllel', 'f8'), ('slope', 'f8'), ('indvMean', 'f8'),
                   ('indvSTD', 'f8'), ('meanFitt', 'f8'), ('stdFitt', 'f8'),
                   ('cvFitMean', 'f8'), ('cvFitSTD', 'f8')])


def plotBoxesGeneMeans(dataArr):
    """ """
    if dataArr.dtype == datype:
        lebls = np.unique(dataArr["Mode"])
    else:
        print("ERROR in plotHistograms(): wrong numpy data type. It should",
              "be:", datype)
        return None
    lbls = []
    for itm in lebls:
        lbls.append(itm.decode())
    ll1 = []
    ll2 = []
    fs = 16
    tkfs = 14
    plt.figure(figsize=(9, 6))
    for itm in lebls:
        ll1.append(dataArr[dataArr["Mode"] == itm]['meanAllel'])
    for itm in lebls:
        ll2.append(dataArr[dataArr["Mode"] == itm]['indvMean'])
    boxprops = dict(linestyle='-', linewidth=2.5, color='k')
    medianprops = dict(linestyle='-', linewidth=2.5)
    whiskerprops = dict(linewidth=2.5)
    capprops = dict(linewidth=2.5)
    flierprops = dict(markersize=10)
    plt.figure(1, figsize=(10, 8))
    plt.boxplot(ll1, labels=lbls, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops)
    plt.ylabel("number of types MHCs in population", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.grid(axis='y')
    plt.figure(2, figsize=(10, 8))
    plt.boxplot(ll2, labels=lbls, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops)
    plt.ylabel("number of MHC in individual", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.grid(axis='y')
    plt.show()


def main():
    """"Main function - the script's main body."""
    if len(sys.argv) <= 1:
        print("One argument is needed:")
        print("  1. Give the name of the data file.")
        sys.exit()
    try:
        dat = np.genfromtxt(sys.argv[1], dtype=datype)
    except:
        print("Cannot load the data from data file. Failed.")
        sys.exit()
    try:
        plotBoxesGeneMeans(dat)
    except:
        print("Failed to plot the data.")
        sys.exit()


if __name__ == "__main__":
    main()
