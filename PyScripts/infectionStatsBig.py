#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the immunocompetence of new (freshly mutated) MHC genes.

Run it e.g. like this:
$ [path_to_scripts]/infectionStatsBig.py [path_to_data_dir] 1500 10 750 'mean'

Running with no arguments gives you tips what they actually should be.

# === side comments ===

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

From Jacek's e-mail from Nov 29 2017:
immunokomp = lm  / (l1p1 + l2p2 + lxpx) gdzie lm - srednia liczba patogenów
zapresentowanych przez nowy allel, lx średnie liczby presentowane przez kolejne
rezydentne allele, px - proporcje rezydentnych alleli w populacji (lepiej
wyskalowane do 100%, czyli nie biorąc pod uwage nowego mutanta)
I tu najlepiej wziąć w momencie pojawienia się mutanta, czyli w pokoleniu 1,
ale wszyskich, które się pojawią (a nie tylko tych, co przetrwają 10 pokoleń)
- więc powstanie boxplot, jak dla prawdopodobieństw rekrutacji, zamiast wykresu
śledzącego w kolejnych pokoleniach

# === end ===

Created on Wed May 10 00:27:47 2017

@author: Piotr Bentkowski :: bentkowski.piotr@gmail.com
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
    failedAttempts = 0
    numbOfFailes = numbOfRand
    if numbOfFailes < 100:
        numbOfFailes = 100
    while(len(geneList) < numbOfRand and failedAttempts <= numbOfFailes):
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
            failedAttempts += 1
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
    """Fed with the MHC ID and the path to InfectionGeneID.csv file it pops a
    two-element Numpy array with the line number and the position in the line
    of the MHC gene in question."""
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
    """Extracts the data from a single, user-specified line. Writes'em to a
    Numpy array."""
    print(firstOccur)
    print(ln.getline(FILE, firstOccur[0][0]))
    ll = re.split(" ", ln.getline(FILE, firstOccur[0][0]))
    del ll[firstOccur[0][1]]
    return np.array(ll[1::], dtype=int)


def totalPlot(geneStats):
    """Does a simple plot to check generic properties of the MHC gene in its
    life."""
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


def statsMHC(mhcID, hostPopSize, avergWay='median', path=os.getcwd()):
    """This is the core function of this script calculating immunocompetence of
    a 'new' MHC gene. It is fed with the ID of an MHC gene and the size of the
    host population. It obtains the numbers generated by the model for that
    particular MHC like the number of hosts that carry it, numbers of presented
    pathogens. Plus the background statistics for the 'resident alleles'.
    Generates a Numpy array containing statistics regarding immunocompetence
    through the genes lifetime."""
    someGene = getIndexGivenGeneID(path + "/InfectionGeneID.csv", mhcID)
    nubPatho = np.array(getDataOfGivenGeneByID(path +
                                               "/InfectionGeneNumbPatho.csv",
                                               someGene), dtype=int)
    hostCount = np.array(getDataOfGivenGeneByID(path +
                                                "/InfectionGeneInHosts.csv",
                                                someGene), dtype=int)
    bkgPathoLd = getBkgroundDataOfGivenGeneByID(path +
                                                "/InfectionGeneNumbPatho.csv",
                                                someGene)
    bkgHostCount = getBkgroundDataOfGivenGeneByID(path +
                                                  "/InfectionGeneInHosts.csv",
                                                  someGene)
    bkgHostWeigth = []
    for itmBkg in bkgHostCount:
        bkgHostWeigth.append(itmBkg / hostPopSize)
    bkg_fitt = []
    for itm in zip(bkgPathoLd, bkgHostWeigth):
        bkg_fitt.append(itm[0] * (itm[1] / hostPopSize))
    meanFitt = np.zeros(len(bkg_fitt))
    if avergWay == 'median':
        for i, sett in enumerate(bkgPathoLd):
            meanFitt[i] = np.median(sett)
    elif avergWay == 'mean':
        for i, sett in enumerate(bkgPathoLd):
            meanFitt[i] = np.mean(sett)
    else:
        print("Select how you want to average the result: 'mean' or 'median'?")
        return None
    meanFitt[meanFitt == 0] = 1.
    return (nubPatho * (hostCount / hostPopSize)) / meanFitt
#    return nubPatho / hostCount


def statsMHC_wrong(mhcID, hostPopSize, avergWay='median', path=os.getcwd()):
    """This is the core function of this script calculating immunocompetence of
    a 'new' MHC gene. It is fed with the ID of an MHC gene and the size of the
    host population. It obtains the numbers generated by the model for that
    particular MHC like the number of hosts that carry it, numbers of presented
    pathogens. Plus the background statistics for the 'resident alleles'.
    Generates a Numpy array containing statistics regarding immunocompetence
    through the genes lifetime. This is the older, probably wrong version."""
    someGene = getIndexGivenGeneID(path + "/InfectionGeneID.csv", mhcID)
    nubPatho = np.array(getDataOfGivenGeneByID(path +
                                               "/InfectionGeneNumbPatho.csv",
                                               someGene), dtype=int)
    hostCount = np.array(getDataOfGivenGeneByID(path +
                                                "/InfectionGeneInHosts.csv",
                                                someGene), dtype=int)
    bkgPathoLd = getBkgroundDataOfGivenGeneByID(path +
                                                "/InfectionGeneNumbPatho.csv",
                                                someGene)
    bkgHostCount = getBkgroundDataOfGivenGeneByID(path +
                                                  "/InfectionGeneInHosts.csv",
                                                  someGene)
    hostWeigth = []
    for itmHo in hostCount:
        hostWeigth.append(itmHo / hostPopSize)
    hostWeigth = np.array(hostWeigth)
    hostWeigth[hostWeigth == 0] = np.nan
    bkgHostWeigth = []
    for itmBkg in bkgHostCount:
        bkgHostWeigth.append(itmBkg / hostPopSize)
    bkg_fitt = []
    for itm in zip(bkgPathoLd, bkgHostWeigth):
        bkg_fitt.append(itm[0] / itm[1])
    meanFitt = np.zeros(len(bkg_fitt))
    if avergWay == 'median':
        for i, sett in enumerate(bkg_fitt):
            meanFitt[i] = np.median(sett)
    elif avergWay == 'mean':
        for i, sett in enumerate(bkg_fitt):
            meanFitt[i] = np.mean(sett)
    else:
        print("Select how you want to average the result: 'mean' or 'median'?")
        return None
    meanFitt[meanFitt == 0] = 1.
    return nubPatho / (hostWeigth * meanFitt)
#    return nubPatho / hostCount


def calcRelatFittManyMHC(geneList, hostPopSize, averWay='median',
                         path=os.getcwd()):
    """Runs the statsMHC() analysis for a long list of preselected MHCs."""
    mhcStatList = []
    for mhc in geneList:
        mhcStatList.append(statsMHC(mhc, hostPopSize, averWay, path))
    return mhcStatList


def removeShortLivedMHC(mhcStatList, minGeneAge=1):
    """Filters out MHC genes that lasted in the population forless then
    minGeneAge of host generations."""
    # Filter out the ones that lasted too few generations
    totalGeneCount = float(len(mhcStatList))
    mhcStatList_trimm = []
    for itm in mhcStatList:
        maxX = len(itm)
        if(maxX >= minGeneAge):
            mhcStatList_trimm.append(itm)
    genesRetained = len(mhcStatList_trimm)
    fracGenesRet = float(genesRetained) / totalGeneCount
    print("There are", genesRetained, "MHC genes that lasted at",
          "least", minGeneAge, "host generations. It is %.2f" % fracGenesRet,
          "of genes submitted for the analysis.")
    return mhcStatList_trimm, totalGeneCount, float(genesRetained)


def getTheMeanRelatFitt(mhcStatList_trimm):
    """Flips the data so the first Numpy array in the output list (index 0)
    contains only the data from the first (index 0) generation of ALL the MHC
    genes (yes, their life is aligned to the same scale). It helps to calculate
    the means in the next step. """
    maxaxX = 0
    for itm in mhcStatList_trimm:
        maxX = len(itm)
        if maxX > maxaxX:
            maxaxX = maxX
    statList = []
    for i in range(maxaxX):
        ll = []
        for itm in mhcStatList_trimm:
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
        plt.plot(geneStats[0][3]/geneStats[0][2],
                 color=(0.75, 0.75, 0.75, 0.75))
    plt.grid(True)
    plt.xlabel("age of a new MHC variant\n[number of host generation from"
               + " de novo mutation]")
    plt.ylabel("number of pathogens presented / number\nof host with the"
               + " new mutation")
    plt.show()


def calcAverageForOneRun(mhcStatList, avergWay='median'):
    """Calculates the median immunocompetence from a single run. Uses the
    precalculated list of individual MCH immunocompetence."""
    # flip the data to time-orientated list
    meanStats = getTheMeanRelatFitt(mhcStatList)
    avrg = np.zeros(len(meanStats))
    if avergWay == 'median':
        for i, itm in enumerate(meanStats):
            avrg[i] = np.median(itm[~np.isnan(itm)])
    elif avergWay == 'mean':
        for i, itm in enumerate(meanStats):
            avrg[i] = np.mean(itm[~np.isnan(itm)])
    else:
        print("Select how you want to average the result: 'mean' or 'median'?")
        return None
#    np.savetxt("avrg.txt", avrg)
    return avrg


def semLogPlotManyMhcStat(mhcStatList, minMhcAge=1, avgWay='median',
                          maxxy=(100, 10e5)):
    """This is the main function that does the plotting. The blue line is the
    median immunocompetence of new MHC introduced by mutations. Uses log-scale
    on the Y axis."""
    fs = 16
    tkfs = 14
    maxaxX = 0
    mhcStatList = removeShortLivedMHC(mhcStatList, minMhcAge)[0]
    ww = calcAverageForOneRun(mhcStatList, avgWay)
    plt.figure(1, figsize=(12, 8))
    for itm in mhcStatList:
        maxX = len(itm)
        if maxX > maxaxX:
            maxaxX = maxX
        plt.semilogy(itm + 1., color=(0.75, 0.75, 0.75, 0.75))
    plt.semilogy(np.ones(maxaxX) * 2., 'k--')
    plt.semilogy(ww + 1., 'b-', lw=2)
    plt.grid(True)
    plt.xlabel("time [host generations after mutation apperence]", fontsize=fs)
    plt.ylabel("mutant’s relative immunocompetence, $log(y + 1)$ ",
               fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    if(maxxy[0] > 10):
        plt.xlim((0, maxxy[0]))
        plt.ylim((0, maxxy[1]))
    plt.show()


def plotManyMhcStat(mhcStatList, minMhcAge=1, avgWay='median',
                    maxxy=(100, 10e5)):
    """This is the main function that does the plotting. The blue line is the
    median immunocompetence of new MHC introduced by mutations. Uses linear
    scale on both axis."""
    fs = 16
    tkfs = 14
    maxaxX = 0
    mhcStatList = removeShortLivedMHC(mhcStatList, minMhcAge)[0]
    ww = calcAverageForOneRun(mhcStatList, avgWay)
    plt.figure(1, figsize=(12, 8))
    for itm in mhcStatList:
        maxX = len(itm)
        if maxX > maxaxX:
            maxaxX = maxX
        plt.plot(itm, color=(0.75, 0.75, 0.75, 0.75))
    plt.plot(np.ones(maxaxX), 'k--')
    plt.plot(ww, 'b-', lw=2)
    plt.grid(True)
    plt.xlabel("time [host generations after mutation apperence]", fontsize=fs)
    plt.ylabel("mutant’s relative immunocompetence",
               fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    if(maxxy[0] > 10):
        plt.xlim((0, maxxy[0]))
        plt.ylim((0, maxxy[1]))
    plt.show()


def main():
    """The main body of the script."""
    if len(sys.argv) <= 5:
        print("4 arguments needed:")
        print("1. Give the name of working directory. May be '.'")
        print("2. Host generation after which stats will be obtained")
        print("3. Minimum number of generations an MHC must exist")
        print("4. Number of genes that need to be selected for analysis")
        print("5. Way to average the results: 'mean' or 'median'")
        sys.exit()
    if(os.path.exists(sys.argv[1])):
        path = sys.argv[1] + "/"
    else:
        print("Data directory non-existent or inaccessible")
        sys.exit()
    try:
        geneList = pickRandomMHCs(path + "InfectionGeneID.csv",
                                  int(sys.argv[2]), int(sys.argv[4]))
        paramsFile = os.path.join(path, 'InputParameters.csv')
        ll = re.split(" ", ln.getline(paramsFile, 7))
        hostPopSize = float(ll[2])
    except Exception:
        print("Cannot load the data. Maybe they're gone...")
        sys.exit()
    xMax = 30
#    mhcStatList = calculateRelatFittManyMHC(geneList, hostPopSize, 'mean')

    mhcStatList = calcRelatFittManyMHC(geneList, hostPopSize, sys.argv[5])

    semLogPlotManyMhcStat(mhcStatList, int(sys.argv[3]), sys.argv[5],
                          (xMax, 10e5))

#    semLogPlotManyMhcStat(mhcStatList, int(sys.argv[3]), 'median',
#                          (xMax, 10e5))

#    plotManyMhcStat(mhcStatList, int(sys.argv[3]), sys.argv[5], (xMax, 200))


if __name__ == "__main__":
    main()
