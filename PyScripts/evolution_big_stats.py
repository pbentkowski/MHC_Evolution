#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 02:38:00 2016

@author: piotr
"""
import re
import os
import linecache as ln
import sys
import numpy as np
import matplotlib.pyplot as plt
import bitstring as bts


outType = np.dtype([('VAR', 'f8'), ('VARX', 'f8'), ('timeMean', 'f8'),
                    ('timeMedian', 'f8'), ('sourceDir', 'S99')])


def lookForVAR(template):
    """Checks which parameters are designated to be investigated as independent
    variables. Gets their line numbers in the file with parameter
    description."""
    varrs = {"VAR": np.nan, "VARX": np.nan}
    for ii, itm in enumerate(template):
        if itm == "VAR":
            varrs["VAR"] = ii
        elif itm == "VARX":
            varrs["VARX"] = ii
        else:
            pass
    return varrs


def loadParamSettings(filepath):
    """Loads model's parametrisation from i.g. InputParameters.csv file into
    a handy list. """
    try:
        paramzList = []
        with open(filepath, 'r') as f:
            for ii, line in enumerate(f):
                if re.search("#", line) or line == "":
                    pass
                else:
                    try:
                        paramzList.append(line.split()[2])
                    except:
                        pass
        return paramzList
    except:
        print("ERROR in loadParamSettings(): Cannot load params into a list.")
        return None


def compareParams(template, paramz):
    """Compares parameters of two runs. They have to be loaded into a list
    first, i.g. with loadParamSettings() function."""
    same = None
    if len(template) == len(paramz):
        for ii, itm in enumerate(zip(template, paramz)):
            try:
                ITM_0 = float(itm[0])
                ITM_1 = float(itm[1])
            except:
                ITM_0 = str(itm[0])
                ITM_1 = str(itm[1])
            if itm[0] == "VARX" or itm[0] == "VAR" or ii <= 1 or ii >= 19:
                pass
            elif ITM_0 == ITM_1:
                same = True
            else:
                same = False
                break
    else:
        print("ERROR in compareParams(): Params lists have different length.")
        return same
    if same is None:
        print("ERROR in compareParams(): Comparison failed to commence.")
    return same


def loadHostPopulation(FILE):
    '''Takes the file with the Host population HostGenomesFile.XXXX.csv and
    picks unique genes from it. Produces two lists: one containing ancestry of
    each gene (tags of all predecessors) and times when each mutation arose in
    the timeline.'''
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
    """ """
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
    """ """
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
    """ """
    ll = []
    for itm in Mut_tags:
        ll.append(len(itm))
    return np.array(ll)


def findMRCA(Mut_tags, Mut_times):
    """Finds the tag, time stamp and index of the most recent common ancestor
    gene from the list of all genes at the population snapshot."""
    if len(findTheOnesAtBeginning(Mut_tags, 0)) != 1:
        print("The most recent common ancestor cannot be established.",
              "There is more than one ancestral gene at the root.")
        return None, np.nan, np.nan
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
    return theMRCAtag, int(Mut_times[0][ii]), ii


def timeOfExistence(Mut_tags, Mut_times):
    """ """
    times = []
    for itm in Mut_times:
        for ii in range(1, len(itm)):
            times.append(int(itm[ii]) - int(itm[ii-1]))
        return np.array(times)


def findLeaves(LIST):
    """ """
    ll = []
    maxLen = 0
    for itm in LIST:
        ll.append((itm[-1], len(itm)))
        if maxLen < len(itm):
            maxLen = len(itm)
    return ll, maxLen


def transTagsToNumpyArr(tagList):
    """ """
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
    """ """
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
    """ """
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
    """ """
    dd = {}
    for item in tagArr:
        dd[np.max(item)] = maxTime
    return dd


def lastingTimeOfGene(timeArr):
    """ """
    lastings = np.zeros(timeArr.shape[0])
    for i, item in enumerate(timeArr):
        for j, ii in enumerate(item):
            if ii < 0 or ii == item[-1]:
                lastings[i] = item[j-1] - item[j-2]
                break
    return lastings


def plotTheTimes(tagArr, timeArr, maxTime, genePairs, maxTimeGenDict, dirrName,
                 linesWdth=(2, 1, 3, 20)):
    """ """
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


def processDataOneFile(FILE):
    """ """
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
    if findMRCA(mutTags, mutTimes):
        mutTimes.sort(key=len, reverse=True)
        mutTags.sort(key=len, reverse=True)
        npMutTags = transTagsToNumpyArr(mutTags)
        npMutTimes = transTimesToNumpyArr(mutTimes, maxTime)
        genePairs, geneTimez = setPairedOriginTags(npMutTags, npMutTimes)
        lastGeneDict = maxGeneLifeDict(npMutTags, maxTime)
        return npMutTags, npMutTimes, maxTime, genePairs, lastGeneDict, \
            geneTimez
    else:
        print("Sorry... No single MRCA.")
        return None


def serchTheDirs(FILE, template, dirr=os.getcwd()):
    """ """
    vv = lookForVAR(template)
    datOut = []
    dataOrdering = ['VAR', 'VARX', 'timeMean', 'timeMedian']
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if filepath == os.path.join(dirName, FILE):
                try:
                    paramzList = loadParamSettings(os.path.join(dirName,
                                                   "InputParameters.csv"))
                except:
                    print("Cannot load the parameters. in dir", dirName)
                    continue
                if compareParams(template, paramzList):
                    try:
                        DATA = processDataOneFile(filepath)
                    except:
                        print("Cannot load the data. in dir", dirName)
                        continue
                    plotTheTimes(DATA[0], DATA[1], DATA[2], DATA[3], DATA[4],
                                 dirName)
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    geneLasts = lastingTimeOfGene(DATA[1])
                    timeMean = np.mean(geneLasts)
                    timeMedian = np.median(geneLasts)
                    datOut.append((var, varx, timeMean, timeMedian, dirName))
    datOut = np.array(datOut, dtype=outType)
    return np.sort(datOut, order=dataOrdering)


def main():
    """ """
    """Main function - the script's main body."""
    if len(sys.argv) <= 2:
        print("Two arguments are needed:")
        print("  1. Give the path to template file.")
        print("  2. Give the name of the output file.")
        sys.exit()
    headerr = 'VAR VARX timeMean timeMedian sourceDir'
    outputFile = str(sys.argv[2])
    try:
        template = loadParamSettings(sys.argv[1])
    except:
        print("Cannot load the template file. Exiting.")
        sys.exit()
    try:
        theData = serchTheDirs("HostGenomesFile.5000.csv", template)
    except:
        print("Failed to process the data. Some serious issues arose.")
        sys.exit()
    if len(theData):
        FMT = '%.4e %.4e %.4e %.4e %s'
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
