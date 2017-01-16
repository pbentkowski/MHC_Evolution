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
# import linecache as ln
import numpy as np
import matplotlib.pyplot as plt
import bitstring as bts


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


def plotTheTimes(tagArr, timeArr, maxTime, genePairs, maxTimeGenDict,
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
    plt.savefig("hitoryOfGenes.png")
    plt.show()


def processDataOneFile(FILE):
    """ """
    try:
        filepath = os.path.join(os.getcwd(), 'InputParameters.csv')
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
#        plotTheTimes(npMutTags, npMutTimes, maxTime, genePairs, lastGeneDict)
        return npMutTags, npMutTimes, maxTime, genePairs, lastGeneDict
    else:
        print("Sorry... No single MRCA.")
        return None


def serchTheDirs(FILE, dirr=os.getcwd()):
    """ """
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if filepath == os.path.join(dirName, FILE):
                try:
                    DATA = processDataOneFile(FILE)
                except:
                    print("Cannot load the data. in dir", dirName)
                    continue
                plotTheTimes(DATA[0], DATA[1], DATA[2], DATA[3], DATA[4])


def main():
    """ """
    try:
        DATA = processDataOneFile(sys.argv[1])
    except:
        print("Cannot load the data. Script aborted.")
        sys.exit()
    if DATA:
        plotTheTimes(DATA[0], DATA[1], DATA[2], DATA[3], DATA[4])
    else:
        print("ERROR in data. Script aborted.")
        sys.exit()


if __name__ == "__main__":
    main()
