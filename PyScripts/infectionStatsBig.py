#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:27:47 2017

ryc 3 w Ejsmond i Radwan 2015, na poczatek może być tak jak D-E
(immunocompetence), ale jak latwiej recruitment probability, to tez może być.
Nie szkodzi, że więcej genów, chcemy wiedzieć, jaka przewagę mają nowe
warianty w róznych kontekstach (zestawach parametrów).

Fig 3.
Presentation ability and probability of recruitment of a new mutant MHC
allele. (A-C) Data points indicate the probability that a mutant allele will
stay in the population for at least 10 host generations (cut-off points longer
or shorter than 10 generations (e.g., 2 or 20) showed similar patterns). (D-F)
Data points indicate the relative pathogen-recognition ability of a mutant MHC
allele relative to the immunocompetence of resident alleles. The average
immunocompetence of resident alleles was measured as the number of pathogens
recognized in a given generation, weighted by the frequency of each resident
allele. Note that a presentation spectrum equal to 1 indicates the threshold at
which a mutant allele is able to present, on average, the same proportion of
pathogens as resident alleles are. (A-F) Characteristics were calculated across
the last 6000 generations of 10 replicates. Scenario labels: HA–heterozygote
advantage, RQ–Red


From the paper fig. 4:
Panels in the bottom row show the dynamics of median mutant allele
immunocompetence in relation to the immunocompetence of resident alleles,
measured as the number of pathogens presented in a given generation (weighted
by the frequency of resident alleles). The measure was calculated across the
last 6000 generations of 10 replicate simulations with the host population size
N = 5000. Note that a presentation spectrum equal to 1 in the bottom row
panels indicate the threshold at which a mutant allele is able to recognize,
on average, the same proportion of pathogens as the resident alleles. For
presentation purposes the scale on the Y-axis in the bottom row panels was
log-transformed. Scenario labels: HA–heterozygote advantage, RQ–Red Queen,
HA+RQ–both.

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
    """Creates a list of unique MHCs IDs form the InfectionGeneID.csv file.
    User can define the generation from which to start recording unique MHCs.
    """
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
    """Fed with the MHC's geneID and the path to InfectionGeneID.csv it creates
    a list of generation numbers when the gene was in existence and its
    position in  the line. Needed for retrieving other data from simulations.
    """
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
    """Fed a data file and an index array created by the function
    getIndexGivenGeneID() it retrieves the data for a given MHC gene."""
    dataInFile = []
    for itm in idxArr:
        ll = re.split(" ", ln.getline(FILE, itm[0]))
        dataInFile.append((int(ll[itm[1]])))
    return dataInFile


def pickRandomMHCs(FILE, firstGen=2, numbOfRand=100):
    """ Creates a list of randomly selected MHC of a given length. User can
    define the host generation to start with."""
    maxLine = file_len(FILE)
    geneList = []
    while(len(geneList) < numbOfRand):
        try:
            ll = re.split(" ", ln.getline(FILE,
                                          rnd.randrange(firstGen, maxLine)))
            mhc = ll[rnd.randrange(1, len(ll)-1)]
        except Exception:
            print("Failed to pick a gene")
            continue
        if mhc not in geneList:
            geneList.append(mhc)
        else:
            continue
    return geneList


def getMHCsStats(randomMHCs):
    """Fed with a list of MHC geneID it creates a list containing MHCs basic
    stats about presence and abundance across host generations."""
    geneStats = []
    for mhc in randomMHCs:
        idxArr = getIndexGivenGeneID("InfectionGeneID.csv", mhc)
        numPath = np.array(getDataOfGivenGeneByID("InfectionGeneNumbPatho.csv",
                                                  idxArr), dtype=int)
        numGenCopy = np.array(getDataOfGivenGeneByID("InfectionGeneNumb.csv",
                                                     idxArr), dtype=int)
        numbHosts = np.array(getDataOfGivenGeneByID("InfectionGeneInHosts.csv",
                                                    idxArr), dtype=int)
        geneStats.append((mhc, numGenCopy, numbHosts, numPath))
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


def getBkgroundDataOfGivenGeneByID(FILE, idxArr):
    """Fed a data file and an index array created by the function
    getIndexGivenGeneID() it retrieves the background data for a given MHC
    gene, what means stats of all the other genes excluding the queried one."""
    dataInFile = []
    for itm in idxArr:
        ll = re.split(" ", ln.getline(FILE, itm[0]))
        del ll[itm[1]]
        dataInFile.append(np.array(ll[1::], dtype=np.int))
    return dataInFile


def getStatsFromOneLine(FILE, firstOccur):
    """ """
    print(firstOccur)
    print(ln.getline(FILE, firstOccur[0][0]))
    ll = re.split(" ", ln.getline(FILE, firstOccur[0][0]))
    del ll[firstOccur[0][1]]
    return np.array(ll[1::], dtype=int)


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
    bkgPathoLoad = getBkgroundDataOfGivenGeneByID("InfectionGeneNumbPatho.csv",
                                                  someGene)
    bkgHostCount = getBkgroundDataOfGivenGeneByID("InfectionGeneInHosts.csv",
                                                  someGene)
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


def getTheMeanRelatFitt(mhcStatList, minGeneAge=1):
    """ """
    maxaxX = 0
    for itm in mhcStatList:
        maxX = len(itm)
        if(maxX < minGeneAge):
            del itm
#            print(maxX)
        if maxX > maxaxX:
            maxaxX = maxX
    statList = []
    for i in range(maxaxX):
        ll = []
        for itm in mhcStatList:
            try:
                ll.append(itm[i])
            except Exception:
                pass
        statList.append(np.array(ll))
    return statList


def randomMhcPlot(geneList, plottedMHCs=23):
    """Test plost of absolute immunocompetence of a set of randomly selected
    MHCs."""
    plt.figure(1, figsize=(8, 6))
    for i in range(plottedMHCs):
        mhcID = geneList[np.random.randint(len(geneList))]
        geneStats = getMHCsStats([mhcID])
        plt.plot(geneStats[0][3]/geneStats[0][2])
        plt.grid(True)
        plt.xlabel("age of a new MHC variant\n[number of host generation from"
                   + " de novo mutation]")
        plt.ylabel("number of pathogens presented / number\nof host with the"
                   + " new mutation")
        plt.show()


def plotManyMhcStat(mhcStatList, minMhcAge=1, maxxy=(100, 10e5)):
    """ """
    fs = 16
    tkfs = 14
    maxaxX = 0
    meanStats = getTheMeanRelatFitt(mhcStatList, minMhcAge)
    ww = np.zeros(len(meanStats))
    for i, itm in enumerate(meanStats):
        ww[i] = np.median(itm[~np.isnan(itm)])
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
    plt.xlabel("time [host generations after mutation apperence]", fontsize=fs)
    plt.ylabel("mutant’s relative immunocompetence, $log(y + 1)$ ",
               fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    if(maxxy[0] > 100):
        plt.xlim((0, maxxy[0]))
        plt.ylim((1, maxxy[1]))
    plt.show()


def main():
    """ """
    if len(sys.argv) <= 4:
        print("4 arguments needed:")
        print("1. Give the name of working directory. May be '.'")
        print("2. Host generation after which stats will be obtained")
        print("3. Minimum number of generations an MHC must exist")
        print("4. Number of genes that need to be selected for analysis")
        sys.exit()
    if(os.path.exists(sys.argv[1])):
        path = sys.argv[1] + "/"
    else:
        print("Data directory non-existent or inaccessible")
        sys.exit()
    try:
        geneList = pickRandomMHCs(path + "InfectionGeneID.csv",
                                  int(sys.argv[2]), int(sys.argv[4]))
#        print(geneList)
    except Exception:
        print("Cannot load the data. Maybe they're gone...")
        sys.exit()
    mhcStatList = calculateRelatFittManyMHC(geneList)
    plotManyMhcStat(mhcStatList, int(sys.argv[3]), (300, 10e5))


if __name__ == "__main__":
    main()
