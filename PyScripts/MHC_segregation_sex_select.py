#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Calculates if MHC genes are correlated in co-occurrence in individuals as
a result e.g. sex selection

Created on Fri Aug 19 18:22:31 2016
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import re
import numpy as np
import bitstring as bts
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib


def loadHostPopulation(FILE):
    '''Takes the file with all the hosts data loads it to a list. Each
    individual is loaded as a list of bit strings.And the population is a list
    of individuals.'''
    LL = []
    ll = []
    nextHost = False
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    if nextHost:
                        LL.append(ll)
                        ll = []
                else:
                    nextHost = True
                    bitt = line.split()[0]
                    ll.append(bts.BitString(bin=bitt).int)
        LL.append(np.array(ll, dtype=int))
#        print j + 1
        return LL
    except IOError as e:
        print("I/O error({0}) in".format(e.errno),
              "loadHostPopulation(): {0}".format(e.strerror))


def calculateCooccurrMatrix(hostList):
    """Takes the list of hosts and their MHC genes and calculates co-occurrence
    of genes found in the population."""
    if len(hostList) == 0:
        print("ERROR in calculateCooccurrMatrix(): the hosts list is empty.")
        return None
    # === fetch unique MHC genes ===
    geneList = []
    for indv in hostList:
        for gene in indv:
            if gene in geneList:
                pass
            else:
                geneList.append(gene)
    geneVec = np.array(geneList, dtype=int)
    N = geneVec.shape[0]
    # === fetch frequencies of occurrence of MHC genes ===
    geneFreq = np.zeros(N)
    for k, gg in enumerate(geneVec):
        cc = 0
        for indv in hostList:
            if(gg in indv):
                cc += 1
        geneFreq[k] = cc
    # === calculate MHC genes co-occurrence matrix ===
    geneFreqMtx = np.zeros((N, N))
    cooccMtx = np.zeros((N, N))
    for i, one in enumerate(geneVec):
        for j, two in enumerate(geneVec):
            count = 0
            for indv in hostList:
                if(one in indv and two in indv):
                    count += 1
            cooccMtx[i, j] = count
            geneFreqMtx[i, j] = np.max([geneFreq[i], geneFreq[j]])
    # === return normalized co-occurrence matrix ===
    return geneFreq, cooccMtx / geneFreqMtx


def clusterAndPlotMHCgenes(mtx, low=0, high=1):
    """Does the clustering of co-occurring MHC genes using SciPy clustering
    package and plots dendrogram and co-occurrence matrix. Takes on input 
    co-occurrence matrix calculated by calculateCooccurrMatrix()."""
    FontSize = 15
#    TiitlFontSize = 12
    TickSize = 14
    Hpad = 0.2
#    mtx = 1.0 - mtx
    print(mtx)
    yy = sch.linkage(mtx, method='centroid')
    plt.clf()
    plt.axes().set_aspect(4.0)
    plt.figure(1, figsize=(16, 8), frameon=False, dpi=100)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    plt.subplot(gs[0])
    matplotlib.rcParams['lines.linewidth'] = 3.7
    Z = sch.dendrogram(yy, orientation='left')
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
    
    

def plotFrequencies(geneFreq):
    '''Plots in how many individuals each MHC gene is present. Takes MHC gene
    frequency vector on input.'''
    N = geneFreq.shape[0]
    FontSize = 12
#    TitlFontSize = 12
#    TickSize = 14
    plt.figure(2, figsize=(10, 5), frameon=False, dpi=100)
    plt.bar(np.arange(N)-0.5, geneFreq)
    plt.xlim((0, N))
    plt.grid(True)
    plt.xlabel("MHC gene tag", fontsize=FontSize)
    plt.ylabel("Number of host individuals with the gene", fontsize=FontSize)
    

def main():
    """ """
    if(len(sys.argv) <= 1):
        print("Give a path to the file with host population data")
        sys.exit()
    try:
        hostList = loadHostPopulation(sys.argv[1])
    except:
        print("ERROR when loading data. Check if the input file has the right",
              "format.")
        sys.exit()
    try:
        geneFreq, cooccMtx = calculateCooccurrMatrix(hostList)
        clusterAndPlotMHCgenes(cooccMtx)
        plotFrequencies(geneFreq)
        print("=================================")
        print("There are", geneFreq.shape[0], "MHC types in population.")
        print("There are", len(hostList[0]), "loci per individual.")
        print("Co-occurrence matrix's min and max:", np.min(cooccMtx), ":",
              np.max(cooccMtx))
    except:
        print("ERROR when processing data. Probably clustering went wrong.",
              "Sorry.")
        sys.exit()
    plt.show() 


if __name__ == "__main__":
    main()
