#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Walks the directory tree looking for model runs which are characterized by same
parametrisation as the template file provided by the user. Then process these
results by fancy stats and plots that processed output on a nice graph.

Created on Tue May 10 18:53:47 2016
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import os
import re
import sys
import datetime as dt
import linecache as ln
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt


# """Data type for loading data from files HostsGeneDivers.csv"""
inType = np.dtype([('time', np.int), ('pop_size', np.int),
                   ('tot_num_of_genes', np.int), ('num_of_MHC_types', np.int),
                   ('Shannon_indx', np.float), ('mean_fitness', np.float),
                   ('std_fitness', np.float)])
# """Data type for storing processed data"""
outType = np.dtype([('VAR', 'f8'), ('VARX', 'f8'), ('meanAllel', 'f8'),
                    ('stdAllel', 'f8'), ('slope', 'f8'), ('indvMean', 'f8'),
                    ('indvSTD', 'f8'), ('meanFitt', 'f8'), ('stdFitt', 'f8'),
                    ('cvFitMean', 'f8'), ('cvFitSTD', 'f8'),
                    ('sourceDir', 'S99')])


def getVarxLabel(filepath):
    """Acquires the name of the VARX variable from input file to later place it
    as the plot's x-axis label."""
    try:
        with open(filepath, 'r') as f:
            for ii, line in enumerate(f):
                if re.search("VARX", line):
                    return line.split()[0]
                else:
                    pass
        print("ERROR in getVarxLabel(): Cannot find the X label.")
        return None
    except:
        print("ERROR in getVarxLabel(): Cannot get the X label.")
        return None


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


def lookForVAR(template):
    """Checks which parameters are designated to be investigated as independent
    variables. Gets their line numbers in the file with parameter
    description."""
    varrs = {"VAR": 0, "VARX": 0}
    for ii, itm in enumerate(template):
        if itm == "VAR":
            varrs["VAR"] = ii
        elif itm == "VARX":
            varrs["VARX"] = ii
        else:
            pass
    return varrs


def readDate(string):
    """Takes a string and tries to convert it into a date. String has to have
    the ISO yyyy-mm-dd format"""
    try:
        dd = string.split("-")
        theDate = dt.date(int(dd[0]), int(dd[1]), int(dd[2]))
        return theDate
    except:
        print("ERROR in readDate(): Bad string format! It has to be ISO's",
              "yyyy-mm-dd format!")
        return None


def loadTheDateFromParamFile(filePar):
    """Takes InputParameters.csv file and tries to figure out what day the run
    was started (line 2 in the file)."""
    try:
        l = re.split(" ", ln.getline(filePar, 2))[2].split(".")[0].split("-")
    except:
        print("ERROR in loadTheDate(): Cannot load the Params file. Check if",
              "the path to the params file is correct as well as its name.")
        return None
    try:
        theDay = dt.date(int(l[0]), int(l[1]), int(l[2]))
        return theDay
    except:
        print("ERROR in loadTheDate(): Cannot covert data into the date",
              "format. Check if the data file has the right flavour.")
        return None


def getTheData(theStartDate, template, EqPt=1000, dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type."""
    vv = lookForVAR(template)
    datOut = []
    dataOrdering = ['VAR', 'VARX', 'meanAllel', 'stdAllel', 'slope']
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.csv') and
               loadTheDateFromParamFile(filepath) >= theStartDate):
                paramzList = loadParamSettings(filepath)
                if compareParams(template, paramzList):
                    l = re.split(" ", ln.getline(filepath, 9))
                    path_spp = float(l[2].split()[0])
                    l = re.split(" ", ln.getline(filepath, 12))
                    pathoNorm = float(l[2].split()[0]) * path_spp
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    dataFilePath = os.path.join(dirName, "HostsGeneDivers.csv")
                    data = np.genfromtxt(dataFilePath, dtype=inType)
                    c0, c1 = poly.polyfit(data['time'][EqPt::],
                                          data['num_of_MHC_types'][EqPt::], 1)
                    meanAlle = data['num_of_MHC_types'][EqPt::].mean()
                    stdAlle = data['num_of_MHC_types'][EqPt::].std()
                    meanFitt = data['mean_fitness'][EqPt::].mean() / pathoNorm
                    stdFitt = np.std(data['mean_fitness'][EqPt::] / pathoNorm)
                    cvFitt = data['std_fitness']/data['mean_fitness']
                    cvFittMean = np.mean(cvFitt[EqPt::]) / pathoNorm
                    cvFittSTD = np.std(cvFitt[EqPt::]) / pathoNorm
                    dataFilePath = os.path.join(dirName,
                                                "HostMHCsNumbUniq_ChrOne.csv")
                    hgsUNIQ = np.genfromtxt(dataFilePath)
                    indvMean = np.mean(hgsUNIQ[EqPt:, 1:])
                    indvSTD = np.std(hgsUNIQ[EqPt:, 1:])
                    datOut.append((var, varx, meanAlle, stdAlle, c1,
                                   indvMean, indvSTD, meanFitt, stdFitt,
                                   cvFittMean, cvFittSTD, dirName))
    datOut = np.array(datOut, dtype=outType)
    return np.sort(datOut, order=dataOrdering)


def buildStats(theData):
    """Averages the data generated by getTheData() function and formats them
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
        meanAll = np.mean(ww['meanAllel'])
        stdAll = np.sqrt(np.sum(ww['stdAllel']**2) / float(len(ww)))
        meanIndv = np.mean(ww['indvMean'])
        stdIndv = np.sqrt(np.sum(ww['indvSTD']**2) / float(len(ww)))
        meanFitt = np.mean(ww['meanFitt'])
        stdFitt = np.sqrt(np.sum(ww['stdFitt']**2) / float(len(ww)))
        meanCvFit = np.mean(ww['cvFitMean'])
        stdCvFit = np.sqrt(np.sum(ww['cvFitSTD']**2) / float(len(ww)))
        meanResult.append((ii[0], ii[1], meanAll, stdAll, meanIndv, stdIndv,
                           meanFitt, stdFitt, meanCvFit, stdCvFit))
    return np.array(meanResult)


def plotAllAllesInPop(meanResult, x_label, logsc='linear'):
    """Uses the array generated by buildStats() function and plots a fancy plot
    of averaged data."""
    FS = 18
    annoSize = int(0.85*FS)
    ll = []
    maxX = 1.15 * float(np.max(meanResult[:, 1]))
    limitz = (0., maxX)
    figSize = (10, 7)
    for itm in meanResult:
        if itm[0] in ll:
            pass
        else:
            ll.append(itm[0])
    # First plot - unique MHC alleles in population
    plt.figure(1, figsize=figSize)
    for var in ll:
        ww = meanResult[meanResult[:, 0] == var]
        plt.errorbar(ww[:, 1], ww[:, 2], ww[:, 3], lw=2, marker="o", ms=8)
        plt.annotate(str(var), xy=(ww[-1, 1], ww[-1, 2]), size=annoSize)
    plt.xlabel(str(x_label), fontsize=FS)
    plt.ylabel("mean number of unique MHC\nalleles in population",
               fontsize=FS)
#    plt.xlim(limitz)
    plt.ylim(ymin=0)
    plt.xscale(logsc)
    plt.tick_params(axis='both', labelsize=annoSize)
    plt.grid(True)
    # Second plot - unique MHC alleles in one chromosome
    plt.figure(2, figsize=figSize)
    for var in ll:
        ww = meanResult[meanResult[:, 0] == var]
        plt.errorbar(ww[:, 1], ww[:, 4], ww[:, 5], lw=2, marker="o", ms=8)
        plt.annotate(str(var), xy=(ww[-1, 1], ww[-1, 4]), size=annoSize)
    plt.xlabel(str(x_label), fontsize=FS)
    plt.ylabel("average number of unique MHC\nalleles in one chromosome",
               fontsize=FS)
#    plt.xlim(limitz)
    plt.ylim(ymin=0)
    plt.xscale(logsc)
    plt.tick_params(axis='both', labelsize=annoSize)
    plt.grid(True)
    plt.figure(3, figsize=figSize)
    for var in ll:
        ww = meanResult[meanResult[:, 0] == var]
        plt.errorbar(ww[:, 1], ww[:, 6], ww[:, 7], lw=2, marker="o", ms=8)
        plt.annotate(str(var), xy=(ww[-1, 1], ww[-1, 6]), size=annoSize)
    plt.xlabel(str(x_label), fontsize=FS)
    plt.ylabel("hosts average fitness normalized per number\nof pathogen" +
               " spp. and pathogen generations", fontsize=FS)
#    plt.xlim(limitz)
    plt.ylim(ymin=0)
    plt.xscale(logsc)
    plt.tick_params(axis='both', labelsize=annoSize)
    plt.grid(True)
    plt.figure(4, figsize=figSize)
    for var in ll:
        ww = meanResult[meanResult[:, 0] == var]
        plt.errorbar(ww[:, 1], ww[:, 8], ww[:, 9], lw=2, marker="o", ms=8)
        plt.annotate(str(var), xy=(ww[-1, 1], ww[-1, 8]), size=annoSize)
    plt.xlabel(str(x_label), fontsize=FS)
    plt.ylabel("hosts average CV fitness normalized per\nnumber of " +
               "pathogen spp. and pathogen generations", fontsize=FS)
#    plt.xlim(limitz)
    plt.ylim(ymin=0)
    plt.xscale(logsc)
    plt.tick_params(axis='both', labelsize=annoSize)
    plt.grid(True)
#    plt.show()


def plotDotMeans(theData):
    """Plots number of MHC alleles in population vs average number of MHC in
    one chromosome."""
    clrs = ['bo', 'go', 'ro', 'co', 'mo', 'yo', 'ko', 'wo']
    clrs += ['bv', 'gv', 'rv', 'cv', 'mv', 'yv', 'kv', 'wv']
    clrs += ['bo', 'go', 'ro', 'co', 'mo', 'yo', 'ko', 'wo']
    FS = 18
    annoSize = int(0.85*FS)
    ll = []
    for itm in theData:
        if (itm['VAR'], itm['VARX']) in ll:
            pass
        else:
            ll.append((itm['VAR'], itm['VARX']))
    plt.figure(5, figsize=(10, 7))
    k = 0
    for var in ll:
        ww = theData[theData['VAR'] == var[0]]
        ww = ww[ww['VARX'] == var[1]]
        lbl = str(var[0]) + " ; " + str(var[1])
        plt.plot(ww['meanAllel'], ww['indvMean'], clrs[k], ms=8, label=lbl)
        k += 1
    print("There are", k, "sets of values")
    plt.legend(loc='lower right', numpoints=1, ncol=2, fontsize=10)
    plt.xlabel("mean number of unique MHC alleles in population",
               fontsize=FS)
    plt.ylabel("average number of unique MHC\nalleles in one chromosome",
               fontsize=FS)
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.tick_params(axis='both', labelsize=annoSize)
    plt.grid(True)
#    plt.show()


def loadTheStuuff(dataSlice, specFile):
    """Loads the data from post-processed files. Useful when working in Ipython
    console.
     . dataSlice - path to DataSlice.csv type of file
     . specFile - path to parameter file type of file"""
    try:
        dd = np.genfromtxt(dataSlice, dtype=outType)
        for itm in dd:
            itm[-1] = itm[-1][2:-1]
        meanResult = buildStats(dd)
        x_Label = getVarxLabel(specFile)
        return dd, meanResult, x_Label
    except:
        print("Bump... Nope...")
        return None


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 2:
        print("Two arguments are needed:")
        print("  1. Give a starting date. It has to be in yyyy-mm-dd format.")
        print("  2. Give the path to template file.")
        sys.exit()
    startDate = None
    headerr = 'VAR VARX meanAllel stdAllel slope indvMean indvSTD meanFitt '\
        + 'stdFitt sourceDir'
    try:
        startDate = readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert argument #1 to a date format.")
        sys.exit()
    if startDate:
        try:
            template = loadParamSettings(sys.argv[2])
            x_Label = getVarxLabel(sys.argv[2])
        except:
            print("Cannot load the template file. Exiting.")
            sys.exit()
        try:
            theData = getTheData(startDate, template)
        except:
            print("Failed to process the data. Some serious issues arose.")
            sys.exit()
        if len(theData):
            FMT = '%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %s'
            open("DataSlice.csv", 'w').close()
            np.savetxt("DataSlice.csv", theData, fmt=FMT, header=headerr,
                       comments='#')
            for itm in theData:
                for ii in range(len(itm) - 1):
                    print(itm[ii], "\t", end=" ")
                print()
            print("Check the output file:", str(os.getcwd()) +
                  "/DataSlice.csv for details.")
            meanResult = buildStats(theData)
            plotAllAllesInPop(meanResult, x_Label)
            plotDotMeans(theData)
            plt.show()
        else:
            print("No data files matching the criterions were found.",
                  "Specify your template file.")
            sys.exit()
    else:
        print("Wrong date format.")
        sys.exit()


if __name__ == "__main__":
    main()
