#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:27:47 2017

@author: piotr
"""
import re
import os
import sys
import linecache as ln
import numpy as np
import random as rnd
import matplotlib.pyplot as plt


def file_len(fname):
    """Get the number of lines in a text file."""
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def getUniqueGenes(FILE, firstIndex=1):
    """ """
    uniqueGenes = []
    with open(FILE) as infile:
        for i, line in enumerate(infile):
            if re.search(r"#", line):
                continue
            elif(i < firstIndex):
                continue
            else:
                LL = line.split()
                for j in range(1, len(LL)):
                    geneID = int(LL[j])
                    if(geneID in uniqueGenes):
                        pass
                    else:
                        uniqueGenes.append(geneID)
    return uniqueGenes


def getIndexGivenGeneID(FILE, geneID, firstIndex=1):
    """ """
    indexInFile = []
    strID = str(geneID)
    with open(FILE) as infile:
        for i, line in enumerate(infile):
            if re.search(r"#", line):
                continue
            elif(i < firstIndex):
                continue
            else:
                if re.search(" " + strID, line):
                    LL = line.split()
                    indexInFile.append((i + 1, LL.index(strID)))
    return np.array(indexInFile, dtype=int)


def getDataOfGivenGeneByID(FILE, idxArr, firstIndex=1):
    """ """
    dataInFile = []
    for itm in idxArr:
        l = re.split(" ", ln.getline(FILE, itm[0]))
        dataInFile.append((int(l[itm[1]])))
    return dataInFile


def pickRandomMHCs(FILE, firstGen=2, numbOfRand=100):
    """ """
    maxLine = file_len(FILE)
    geneList = []
    while(len(geneList) < numbOfRand):
        try:
            l = re.split(" ", ln.getline(FILE,
                                         rnd.randrange(firstGen, maxLine)))
            mhc = l[rnd.randrange(1, len(l)-1)]
        except:
            print("Failed to pick a gene")
            continue
        if mhc not in geneList:
            geneList.append(mhc)
        else:
            continue
    return geneList


def getMHCsStats(randomMHCs):
    """ """
    geneStats = []
    for mhc in randomMHCs:
        idxArr = getIndexGivenGeneID("InfectionGeneID.csv", mhc)
        present = np.array(getDataOfGivenGeneByID("InfectionGeneNumbPatho.csv",
                                                  idxArr), dtype=int)
        abundan = np.array(getDataOfGivenGeneByID("InfectionGeneNumb.csv",
                                                  idxArr), dtype=int)
        geneStats.append((mhc, present, abundan))
    return geneStats


def getFirstOccurence(FILE, geneID, firstIndex=1):
    """ """
    indexInFile = []
    strID = str(geneID)
    with open(FILE) as infile:
        for i, line in enumerate(infile):
            if re.search(r"#", line):
                continue
            elif(i < firstIndex):
                continue
            else:
                if re.search(" " + strID, line):
                    LL = line.split()
                    indexInFile.append((i + 1, LL.index(strID)))
                    break
    if(len(indexInFile) >= 0):
        return np.array(indexInFile, dtype=int)
    else:
        print("MHC ID not found")
        return None


def getTheLinesByIndexes(FILE, idxArr):
    """ """
    dataInFile = []
    for itm in idxArr:
        l = re.split(" ", ln.getline(FILE, itm[0]))
        del l[itm[1]]
        dataInFile.append(np.array(l[1::], dtype=np.int))
    return dataInFile


def getStatsFromOneLine(FILE, firstOccur):
    """ """
    print(firstOccur)
    print(ln.getline(FILE, firstOccur[0][0]))
    l = re.split(" ", ln.getline(FILE, firstOccur[0][0]))
    del l[firstOccur[0][1]]
    return np.array(l[1::], dtype=int)


def totalPlot(geneStats):
    """ """
    for oneMHC in geneStats:
        plt.figure(1)
        plt.subplot(211)
        plt.plot(oneMHC[2])
        plt.ylabel("number of host with\ngiven MHC variant")
        plt.grid(True)
        plt.subplot(212)
        plt.plot(oneMHC[1])
        plt.xlabel("relative time [host generations since a gene emerged]")
        plt.ylabel("number of pathogen presentations\nper host generation")
        plt.grid(True)

        plt.figure(2)
        plt.plot(oneMHC[1]/oneMHC[2])
        plt.xlabel("relative time [host generations since a gene emerged]")
        plt.ylabel("MHC presentation efficiency\n(number of pat." +
                   " presentations per one host")
        plt.grid(True)
    plt.show()


def statsMHC(mhcID, path=os.getcwd()):
    """ """
    someGene = getIndexGivenGeneID("InfectionGeneID.csv", mhcID)
    nubPatho = np.array(getDataOfGivenGeneByID("InfectionGeneNumbPatho.csv",
                                               someGene), dtype=int)
    hostCount = np.array(getDataOfGivenGeneByID("InfectionGeneInHosts.csv",
                                                someGene), dtype=int)
    bkgPathoLoad = getTheLinesByIndexes("InfectionGeneNumbPatho.csv", someGene)
    bkgHostCount = getTheLinesByIndexes("InfectionGeneInHosts.csv", someGene)
    bkg_fitt = []
    for itm in zip(bkgPathoLoad, bkgHostCount):
        bkg_fitt.append(itm[0] / itm[1])
    meanFitt = np.zeros(len(bkg_fitt))
    for i, sett in enumerate(bkg_fitt):
        meanFitt[i] = np.mean(sett)
    meanFitt[meanFitt == 0] = np.nan
    return (nubPatho / hostCount) / meanFitt


def calculateRelatFittManyMHC(geneList):
    """ """
    mhcStatList = []
    for mhc in geneList:
        mhcStatList.append(statsMHC(mhc))
    return mhcStatList


def getTheMeanRelatFitt(mhcStatList):
    """ """
    maxaxX = 0
    for itm in mhcStatList:
        maxX = len(itm)
        if maxX > maxaxX:
            maxaxX = maxX
    statList = []
    for i in range(maxaxX):
        ll = []
        for itm in mhcStatList:
            try:
                ll.append(itm[i])
            except:
                pass
        statList.append(np.array(ll))
    return statList


def plotManyMhcStat(mhcStatList, maxx=100, maxy=10e5):
    """ """
    fs = 16
    tkfs = 14
    maxaxX = 0
    meanStats = getTheMeanRelatFitt(mhcStatList)
    ww = np.zeros(len(meanStats))
    for i, itm in enumerate(meanStats):
        ww[i] = np.mean(itm[~np.isnan(itm)])
#        print(ww[i], end=", ")
    plt.figure(1, figsize=(12, 8))
    for itm in mhcStatList:
        maxX = len(itm)
        if maxX > maxaxX:
            maxaxX = maxX
        plt.semilogy(1. + itm, color=(0.75, 0.75, 0.75, 0.75))
    plt.semilogy(np.ones(maxaxX) * 2., 'k--')
    plt.semilogy(1. + ww, 'b-', lw=2)
    plt.grid(True)
    plt.xlabel("Time since gene emerged [host generations]", fontsize=fs)
    plt.ylabel("gene's relative fitness , $log(x + 1)$ ", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    if(maxx > 100):
        plt.xlim((0, maxx))
        plt.ylim((1, 1e5))
    plt.show()


def main():
    """ """
    if len(sys.argv) < 2:
        print("Give the name of working directory. May be '.'")
        sys.exit()
    if(os.path.exists(sys.argv[1])):
        path = sys.argv[1] + "/"
    else:
        print("Data directory non-existent or inaccessible")
        sys.exit()
    try:
        geneList = pickRandomMHCs(path + "InfectionGeneID.csv", 1000, 500)
#        print(geneList)
    except:
        print("Cannot load the data. Maybe they're gone...")
        sys.exit()
    mhcStatList = calculateRelatFittManyMHC(geneList)
    plotManyMhcStat(mhcStatList, 300)


if __name__ == "__main__":
    main()
