#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Iterates through directories and looks for the file with the Host population
snapshot called HostGenomesFile.XXXX.csv and the InputParameters.csv file with
the parameters used it the run. It extracts information about the genes origin
like ancestry tree and MRCA.

Created on Thu Dec  8 02:38:00 2016

@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import os
import linecache as ln
import sys
import numpy as np
import matplotlib.pyplot as plt
import bitstring as bts
import packed_plots_of_MHC_alleles as ppma


outType = np.dtype([('VAR', 'f8'), ('VARX', 'f8'), ('MRCA_time', 'f8'),
                    ('maxMutNumb', 'f8'), ('numOfGenes', 'f8'),
                    ('sourceDir', 'S99')])


def loadHostPopulation(FILE):
    '''Takes the file with the Host population HostGenomesFile.XXXX.csv and
    picks unique genes from it. Produces two lists: one containing ancestry of
    each gene (tags of all predecessors) and the second, corresponding
    containing times when each mutation arose in the genes time line.'''
    B_list = []
    Mut_tags = []
    Mut_times = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    continue
                else:
                    LL = line.split()
                    bb = bts.BitString(bin=LL[0]).int
                    if bb in B_list:
                        pass
                    else:
                        # print(LL[5::2])
                        B_list.append(bb)
                        tagz = LL[5::2]
                        tagz.append(LL[3])
                        Mut_tags.append(tagz)
                        timez = LL[4::2]
                        timez.append(LL[2])
                        Mut_times.append(timez)
        return Mut_tags, Mut_times
    except IOError as e:
        print("I/O error({0}) in".format(e.errno) +
              " loadTheHostPopulation(): {0}".format(e.strerror))


def loadRawBitstrings(FILE):
    """Loads just the raw bit strings to a list from the Host population
    HostGenomesFile.XXXX.csv file. Used only for debugging.
    """
    B_list = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    continue
                else:
                    LL = line.split()
                    bb = LL[0]
                    if bb in B_list:
                        pass
                    else:
                        B_list.append(bb)
        return B_list
    except IOError as e:
        print("I/O error({0}) in".format(e.errno) +
              " loadTheHostPopulation(): {0}".format(e.strerror))


def findTheOnesAtBeginning(Mut_tags, jj=0):
    """Finds the tag of the first gene in population’s known history."""
    ll = []
    for itm in Mut_tags:
        try:
            if itm[jj] in ll:
                pass
            else:
                ll.append(itm[jj])
        except:
            pass
    return ll


def numberOfMutList(Mut_tags):
    """Creates a list of numbers of mutations each gene had in its history"""
    ll = []
    for itm in Mut_tags:
        ll.append(len(itm))
    return np.array(ll)


def findMRCA(Mut_tags, Mut_times):
    """Finds the tag, time stamp and index of the most recent common ancestor
    gene from the list of all genes at the population snapshot.
    Returns: MRCA gene tag, time of MRCA origin, MRCA index in Mut_tags list,
             time of the MRCA
    """
    if len(findTheOnesAtBeginning(Mut_tags, 0)) != 1:
        print("The most recent common ancestor cannot be established.",
              "There is more than one ancestral gene at the root.")
        return None, np.nan, np.nan, np.nan
    mutNumb = numberOfMutList(Mut_tags)
    maxx = np.max(mutNumb)
    theMRCAtag = Mut_tags[0][0]
    ii = 0
    for x in range(maxx):
        if len(findTheOnesAtBeginning(Mut_tags, x)) == 1:
            theMRCAtag = findTheOnesAtBeginning(Mut_tags, x)[0]
            ii = x
        else:
            break
    minn = np.inf
    kk = ii + 1
    for item in Mut_times:
        try:
            indx = int(item[kk])
        except:
            continue
        if indx < minn:
            minn = indx
    return theMRCAtag, int(Mut_times[0][ii]), ii, minn


def transTagsToNumpyArr(tagList):
    """Takes list of tags created by loadHostPopulation() and transforms in into
    a Numpy array where number of columns equals to the number of unique genes
    and the number of rows equals to the number of mutation events in the
    history of the gene that had the most of these events."""
    maxLen = 0
    for itm in tagList:
        if maxLen < len(itm):
            maxLen = len(itm)
    arr = -1 * np.ones((len(tagList), maxLen), dtype='i8')
    for i, itm in enumerate(tagList):
        for j, ii in enumerate(itm):
            arr[i, j] = ii
    return arr


def transTimesToNumpyArr(timesList, finito):
    """Takes list of mutation times created by loadHostPopulation() and
    transforms in into a Numpy array where number of columns equals to the
    number of unique genes and the number of rows equals to the number of
    mutation events in the history of the gene that had the most of these
    events, 'finito' is the maximal number of generations a.k.a. simulation
    time. """
    maxLen = 0
    for itm in timesList:
        if maxLen < len(itm):
            maxLen = len(itm)
    arr = -1 * np.ones((len(timesList), maxLen+1))
    for i, itm in enumerate(timesList):
        for j, ii in enumerate(itm):
            arr[i, j] = ii
        arr[i, j+1] = finito
    return arr


def setPairedOriginTags(tagArr, timeArr):
    """Takes the Numpy-transformed tag and time arrays and creates a data
    structure where there are only unique ancestor-descendant pairs of genes
    and pairs of corresponding times of gene origin."""
    maxLen = 0
    for itm in tagArr:
        if maxLen < len(itm):
            maxLen = len(itm)
    genePairs = []
    geneTimez = []
    for i in range(maxLen-1):
        tag_ll = []
        time_ll = []
        for j, itm in enumerate(tagArr):
            geneTpl = (itm[i], itm[i+1])
            if itm[i+1] != -1 and geneTpl not in tag_ll:
                    tag_ll.append(geneTpl)
                    time_ll.append((timeArr[j][i+1], timeArr[j][i+2]))
        genePairs.append(tag_ll)
        geneTimez.append(time_ll)
    return genePairs, geneTimez


def maxGeneLifeDict(tagArr, maxTime):
    """Creates an utility data structure: a dictionary of the surviving gene
    tags and the number of host generations (a.k.a. simulation time)"""
    dd = {}
    for item in tagArr:
        dd[np.max(item)] = maxTime
    return dd


def plotTheTimes(tagArr, timeArr, maxTime, genePairs, maxTimeGenDict, dirrName,
                 linesWdth=(2, 1, 3, 20)):
    """Plots the history tree of the population at the snapshot. Saves it to
    file. """
    hLineWidth = linesWdth[0]
    vLineWidth = linesWdth[1]
    lastLineWdth = linesWdth[2]
    FS = linesWdth[3]
    plt.figure(1, figsize=(20, 16))
    lline = 0
    checkList = [-1]
    dd = {}
    for i, itm in enumerate(tagArr):
        for j, ii in enumerate(itm):
            if ii not in checkList:
                checkList.append(ii)
                plt.hlines(lline, timeArr[i][j], timeArr[i][j+1],
                           colors='k', lw=hLineWidth)
                dd[ii] = (timeArr[i][j], lline, timeArr[i][j+1])
                lline += 1
            else:
                pass
    for k, item in enumerate(genePairs):
        for l, pair in enumerate(item):
            y1 = dd[pair[0]][1]
            y2 = dd[pair[1]][1]
            x1 = dd[pair[1]][0]
            x2 = dd[pair[0]][2]
            plt.vlines(x1, y1, y2, colors='k', lw=vLineWidth)
            plt.hlines(y1, x1, x2, colors='k', lw=hLineWidth)
    for item in maxTimeGenDict:
        plt.hlines(dd[item][1], dd[item][0], maxTimeGenDict[item],
                   colors='r', lw=lastLineWdth)
    plt.xlim((0, maxTime))
    plt.ylim(ymin=0)
    plt.xlabel("time [hosts generations]", fontsize=FS)
    plt.xticks(size=FS-2)
    plt.yticks([])
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(dirrName + "/hitoryOfGenes.png")
#    plt.show()
    plt.clf()


def processDataOneFile(FILE):
    """A meta-function that takes the host population snapshot file
    HostGenomesFile.XXXX.csv, also loads some parameters for the file
    InputParameters.csv and pre-processes the population to produce more
    useful data structures."""
    try:
        filepath = os.path.join(os.path.dirname(FILE), 'InputParameters.csv')
        print(filepath)
        l = re.split(" ", ln.getline(filepath, 13))
        maxTime = float(l[2].split()[0])
    except:
        print("Can't load the data from InputParameters.csv.",
              "Check if it exist.")
        return None
    try:
        mutTags, mutTimes = loadHostPopulation(FILE)
    except:
        print("Can't load the host population snapshot file.")
        return None
    mrcaTag, mrcaOri, mrcaIdx, mrcaTime = findMRCA(mutTags, mutTimes)
    if mrcaTag:
        mutTimes.sort(key=len, reverse=True)
        mutTags.sort(key=len, reverse=True)
        npMutTags = transTagsToNumpyArr(mutTags)
        npMutTimes = transTimesToNumpyArr(mutTimes, maxTime)
        genePairs, geneTimez = setPairedOriginTags(npMutTags, npMutTimes)
        lastGeneDict = maxGeneLifeDict(npMutTags, maxTime)
        return npMutTags, npMutTimes, maxTime, genePairs, lastGeneDict, \
            geneTimez, mrcaTime
    else:
        print("Sorry... No single MRCA.")
        return None


def serchTheDirs(FILE, template, dirr=os.getcwd()):
    """Walk the directory tree in search of model runs and process each
    simulation individually. Produces some meta-statistics regarding the
    results geathered in Numpy structured array."""
    vv = ppma.lookForVAR(template)
    datOut = []
    dataOrdering = ['VAR', 'VARX', 'timeMean', 'timeMedian']
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if filepath == os.path.join(dirName, FILE):
                try:
                    paramzList = ppma.loadParamSettings(os.path.join(dirName,
                                                        "InputParameters.csv"))
                except:
                    print("Cannot load the parameters. in dir", dirName)
                    continue
                if ppma.compareParams(template, paramzList):
                    try:
                        DATA = processDataOneFile(filepath)
                    except:
                        print("Cannot load the data. in dir", dirName)
                        continue
                    plotTheTimes(DATA[0], DATA[1], DATA[2], DATA[3], DATA[4],
                                 dirName)
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    datOut.append((var, varx, DATA[6], DATA[0].shape[1],
                                   DATA[0].shape[0], dirName))
    datOut = np.array(datOut, dtype=outType)
    return np.sort(datOut, order=dataOrdering)


def buildStats(theData):
    """Averages the data generated by serchTheDirs() function and formats them
    into Numpy array suitable for latter plotting."""
    LL = []
    for itm in theData:
        tt = (itm['VAR'], itm['VARX'])
        if tt in LL:
            pass
        else:
            LL.append(tt)
    meanResult = []
    for ii in LL:
        ww = theData[theData['VAR'] == ii[0]]
        ww = ww[ww['VARX'] == ii[1]]
#        NN = float(len(ww))
        meanAll = np.mean(ww['MRCA_time'])
        medAll = np.median(ww['MRCA_time'])
        meanResult.append((ii[0], ii[1], meanAll, medAll))
    return np.array(meanResult)


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 2:
        print("Two arguments are needed:")
        print("  1. Give the path to template file.")
        print("  2. Give the name of the output file.")
        sys.exit()
    headerr = 'VAR VARX MRCA_time maxMutNumb numOfGenes sourceDir'
    outputFile = str(sys.argv[2])
    try:
        template = ppma.loadParamSettings(sys.argv[1])
    except:
        print("Cannot load the template file. Exiting.")
        sys.exit()
    try:
        theData = serchTheDirs("HostGenomesFile.5000.csv", template)
    except:
        print("Failed to process the data. Some serious issues arose.")
        sys.exit()
    if len(theData):
        FMT = '%.4e %.4e %.4e %.4e %.4e %s'
        open(outputFile, 'w').close()
        np.savetxt(outputFile, theData, fmt=FMT, header=headerr,
                   comments='#')
        for itm in theData:
            for ii in range(len(itm) - 1):
                print(itm[ii], "\t", end=" ")
            print()
        print("Check the output file:", str(os.getcwd()) + "/" + outputFile +
              " for details.")
    else:
        print("No data files matching the criterions were found.",
              "Specify your template file.")
        sys.exit()


if __name__ == "__main__":
    main()
