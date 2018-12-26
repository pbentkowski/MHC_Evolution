#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Your doc string goes here...

Created on Wed Jul 13 18:26:33 2016
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import os
import sys
import json
import numpy as np
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def loadParams(dirname):
    """ """
    try:
        paramsFile = os.path.join(dirname, 'InputParameters.json')
        with open(paramsFile) as f:
            prms = json.load(f)
        popSize = int(prms['pathogen_population_size'])
        sppNumber = int(prms['number_of_pathogen_species'])
        return popSize, sppNumber
    except Exception:
        print("ERROR in loadParams(): cannot find parameters in file.",
              "Check if the file exists.")
        return None, None


def loadPresentedSpecies(filepath):
    """ """
    LL = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if re.search(" ===", line):
                    LL.append(list(map(int, line.split("are:")[1].
                              split(" ===")[0].split())))
                else:
                    pass
        return LL
    except Exception:
        print("ERROR in loadPresentedSpecies(): Cannot load the presented",
              "pathogen species.")
        return None


def calculateCoocurenceMtx(LL, sppNumber, popSize):
    '''Calculates the matrix of co-presentation of different pathogen species
    by simply counting in how many hosts two species of pathogens were
    presented together.

    ==============  ==========================================================
    argument        info
    ==============  ==========================================================
    `LL`            List with presented pathogens by each single host created
                    by :py:func:`pathogen_spp_cooccur.loadPresentedSpecies`
    `sppNumber`     number of pathogen species
    `popSize`       number of host individuals
    ==============  ==========================================================
    '''
    mtx = np.zeros((sppNumber, sppNumber))
    for i in range(sppNumber):
        for j in range(i+1, sppNumber):
            for itm in LL:
                if i in itm and j in itm:
                    mtx[i, j] += 1
                    mtx[j, i] += 1
    return mtx / float(popSize)


def clusterSpecies(mtx, low=0, high=1):
    ''' '''
    FontSize = 15
#    TitlFontSize = 12
    TickSize = 14
    Hpad = 0.2
#    mtx = 1.0 - mtx
    print(mtx)
    yy = sch.linkage(mtx, method='centroid')
    plt.clf()
    plt.axes().set_aspect(4.0)
    plt.figure(1, figsize=(48, 12), frameon=False, dpi=100)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    plt.subplot(gs[0])
    Z = sch.dendrogram(yy, orientation='right')
    plt.subplot(gs[1])
    index = Z['leaves']
    mtx = mtx[index, :]
    mtx = mtx[:, index]
    plt.pcolormesh(mtx, cmap=plt.cm.jet, vmin=low, vmax=high)
#    plt.title(tittle, fontsize=TitlFontSize)
    plt.axis([0, len(mtx), 0, len(mtx)])
    plt.xticks([])
    plt.yticks([])
    cb = plt.colorbar(format=r"%.1f", orientation='vertical')
    cb.ax.tick_params(labelsize=TickSize)
    cb.set_label(r'distance', fontsize=FontSize)
    cb.ax.minorticks_on()
    plt.tight_layout(h_pad=Hpad)
    plt.show()


def main():
    ''' '''
    if len(sys.argv) <= 1:
        print("Specify the file with hosts genetic data. It is usualy named:",
              "HostGenomesFile.XXXX.csv or similar.")
        sys.exit()
    try:
        popSize, sppNumber = loadParams(os.getcwd())
        LL = loadPresentedSpecies(sys.argv[1])
    except Exception:
        print("Cannot load the data. Check if you are in the right folder.")
        sys.exit()
    print("Number of pathogen spp:", sppNumber)
    print("Host population size  :", popSize)
    mtx = calculateCoocurenceMtx(LL, sppNumber, popSize)
    clusterSpecies(mtx)


if __name__ == "__main__":
    main()
