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
import random


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


def createCloneList(hostList, removeConvergence="no"):
    ''' '''
    cloneList = []
    cloneList.append(hostList[0])
    for indv in hostList:
        newOne = False
        for indvC in cloneList:            
            if np.array_equal(indv, indvC):
                newOne = False
                break   
            else:
                newOne = True
        if newOne:            
            cloneList.append(indv)
    if removeConvergence == "no":
        return cloneList
    else:
        sortedCloneList = []
        for indvS in cloneList:
            sortedCloneList.append(np.sort(indvS))
        reducedCloneList = []
        reducedCloneList.append(sortedCloneList[0])
        for indv in sortedCloneList:
            newOne = False
            for indvC in reducedCloneList:            
                if np.array_equal(indv, indvC):
                    newOne = False
                    break   
                else:
                    newOne = True
            if newOne:            
                reducedCloneList.append(indv)
        return reducedCloneList
    

def calculateCloneSimMatrix(cloneList):
    """ """
    N = len(cloneList)
    simMtx = np.zeros((N, N))
    for i, indv in enumerate(cloneList):
        for j, indv2 in enumerate(cloneList):
            count = 0.
            for gene in indv:
                if gene in indv2:
                    count += 1.
            simMtx[i, j] = count
    return simMtx / float(len(cloneList[0]))
    
        


def calculateGeneCooccurrMatrix(hostList):
    """Takes the list of hosts and their MHC genes and calculates co-occurrence
    of genes found in the population."""
    if len(hostList) == 0:
        print("ERROR in calculateGeneCooccurrMatrix(): the hosts list is empty.")
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


def clusterAndPlotSimMtx(mtx, figgNum=1, low=0, high=1):
    """Does the clustering of co-occurring MHC genes using SciPy clustering
    package and plots dendrogram and co-occurrence matrix. Takes on input 
    co-occurrence matrix calculated by calculateGeneCooccurrMatrix()."""
    FontSize = 15
#    TiitlFontSize = 12
    TickSize = 14
    Hpad = 0.2
#    mtx = 1.0 - mtx
    print(mtx)
    yy = sch.linkage(mtx, method='centroid')
    plt.clf()
    plt.axes().set_aspect(4.0)
    plt.figure(figgNum, figsize=(16, 8), frameon=False, dpi=100)
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
    plt.figure(3, figsize=(8, 5), frameon=False, dpi=100)
    plt.bar(np.arange(N)-0.5, geneFreq)
    plt.xlim((0, N))
    plt.grid(True)
    plt.xlabel("MHC gene tag", fontsize=FontSize)
    plt.ylabel("number of host individuals with the gene", fontsize=FontSize)
    
    
def countAndPlotClonalIndiv(hostList, fracOfPop=0.2):
    """ """
    FontSize = 11
    TickSize = 9
    selectIndiv = int(fracOfPop * float(len(hostList)))
    ww =random.sample(hostList, selectIndiv)
    ll = []
    for i, itm in enumerate(ww):
        count = 0
        for j, itm2 in enumerate(ww):
            if (np.array_equal(itm, itm2) and i != j):
                count += 1
        ll.append(count)
    ll = np.array(ll, dtype=float)
    cloneFreq = ll / float(len(ww))
    # === plot it ===
    plt.figure(3, figsize=(8, 5))  #, frameon=False, dpi=100)
    plt.hist(cloneFreq, normed=True)
    plt.grid(True)
    plt.xlabel("fraction of population size", fontsize=FontSize)
    plt.ylabel("frequency of occurrence", fontsize=FontSize)
    plt.xticks(fontsize=TickSize)
    plt.yticks(fontsize=TickSize)


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
        geneFreq, cooccMtx = calculateGeneCooccurrMatrix(hostList)
        cloneList = createCloneList(hostList, "yes")
        cloneSimMtx = calculateCloneSimMatrix(cloneList)
        clusterAndPlotSimMtx(cooccMtx, 1)
        plt.show()
        clusterAndPlotSimMtx(cloneSimMtx, 1)
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
