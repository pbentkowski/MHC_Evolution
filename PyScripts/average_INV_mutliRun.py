#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
     Your doc string here, please...

Created on Sun Jul 21 15:07:56 2019

@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as subp

import packed_plots_of_MHC_alleles as ppma


def awkMeanINV(path='.'):
    """ """
#    inFile = os.path.join(path, "NumberOfMhcAfterMating.csv")
    inFile = os.path.join(path, "NumberOfMhcBeforeMating.csv")
    outFile = os.path.join(path, "MeanInvdMhcNumb.csv")
    awk = """awk 'NR==1 { next }
               { T=0
                  for(N=2; N<=NF; N++) T+=$N;
                  T/=(NF-1)
                  print $1, T }' """ + inFile + " > " + outFile
    subp.check_output(awk, shell=True)


def loadMeanInvdMhcNumb(path='.'):
    """ """
    theFile = os.path.join(path, "MeanInvdMhcNumb.csv")
    try:
        return pd.read_csv(theFile, delimiter=" ", header=None,
                           names=['time', 'meanINV'])
    except Exception:
        print("ERROR in loadMeanInvdMhcNumb(): Failes to load the",
              "MeanInvdMhcNumb.csv file. Check if it exists")
        return None


def getTheData(theStartDate, templateList, dirr=os.getcwd()):
    """ """
    datOut = []
    vv = ppma.lookForVARinList(templateList)
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.json') and
               ppma.loadTheDateFromParamFile(filepath) >= theStartDate):
                paramzList = ppma.loadParamSettings(filepath)
                if ppma.compareParams(templateList, paramzList):
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    awkMeanINV(dirName)
                    meanINV = loadMeanInvdMhcNumb(dirName)
                    datOut.append((var, varx, meanINV))
                    print("Done dir:", dirName)
    return datOut


def aggrDataByRuns(datOut):
    """ """
    runParm = []
    lenMin = int(1e10)
    for itm in datOut:
        if lenMin < len(itm[2]):
            lenMin = len(itm[2])
        if (itm[0], itm[1]) not in runParm:
            runParm.append((itm[0], itm[1]))
    aggrOut = []
    for rP in runParm:
        ww = []
        for run in datOut:
            if (run[0], run[1]) == rP:
                ww.append(run[2]['meanINV'][:lenMin])
        WW = np.transpose(np.array(ww))
        aggrOut.append((rP, np.mean(WW, axis=1), np.std(WW, axis=1)))
    return aggrOut


def aggrDataByRunsCI(datOut):
    """ """
    runParm = []
    lenMin = int(1e10)
    for itm in datOut:
        if lenMin > len(itm[2]):
            lenMin = len(itm[2])
        if (itm[0], itm[1]) not in runParm:
            runParm.append((itm[0], itm[1]))
    aggrOut = []
    for rP in runParm:
        ww = []
        for run in datOut:
            if (run[0], run[1]) == rP:
                ww.append(run[2]['meanINV'][:lenMin])
        WW = np.transpose(np.array(ww))
        ci = np.zeros(len(WW))
        for i, w in enumerate(WW):
            ci[i] = ppma.confidence_interval(w)
        try:
            aggrOut.append((rP, np.mean(WW, axis=1), ci))
        except Exception:
            print("problem in", rP, np.mean(WW, axis=1))
    return aggrOut


def smoothing(inv, NN=500):
    """ """
    return np.convolve(inv, np.ones((NN,))/NN, mode='valid')


def plotAggrOut(aggrOut, figtag=''):
    """ """
    clrs = ['C0', 'C1', 'C2', 'C3']
    plt.figure(1, figsize=(9, 6))
    for i, agRun in enumerate(aggrOut):
        agrRun = smoothing(agRun[1])
        agrSTD = smoothing(agRun[2])
        x = np.arange(0, len(agrRun))
        plt.plot(x, agrRun, lw=2, color=clrs[i], label=str(int(2*agRun[0][1])))
        plt.fill_between(x, agrRun+agrSTD, agrRun-agrSTD, facecolor=clrs[i],
                         alpha=0.5)
    plt.legend(title="INV at init.")
    plt.grid()
    plt.xlabel("time [host generations]")
    plt.ylabel(r"individual number of MHC variants (INV $\pm$ 95%CI)")
    plt.ylim(bottom=0)
    plt.tight_layout()
    figName = "INV_time_trajc" + figtag + ".png"
    plt.savefig(figName)
    plt.show()


def main():
    """ """
    """Main function - the script's main body."""
    if len(sys.argv) <= 3:
        print("Two arguments are needed:")
        print("  1. Give a starting date. It has to be in yyyy-mm-dd format.")
        print("  2. Give the path to template file.")
        print("  3. Give the output figure name's prefix (e.g. the number",
              "of individual number of MHC variants.")
        sys.exit()
    startDate = None
    try:
        startDate = ppma.readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert argument #1 to a date format.")
        sys.exit()
    if startDate:
        try:
            template = ppma.loadParamSettings(sys.argv[2])
            if template is None:
                print("Failed to load the template file. Exiting.",
                      "Check if the path is correct - you may wish to provide",
                      "an absolute path.")
                sys.exit()
            figLabel = ppma.getVarxLabel(sys.argv[2])
        except Exception:
            print("Cannot load the template file. Exiting.")
            sys.exit()
        if True:
            # third argument is very important
            theData = getTheData(startDate, template)
            print(theData)
#        except Exception:
        else:
            print("Failed to process the data. Some serious issues arose.")
            sys.exit()
        if len(theData):
            aggrOutCI = aggrDataByRunsCI(theData)
            plotAggrOut(aggrOutCI, figLabel)
        else:
            print("No data files matching the criterions were found.",
                  "Specify your template file.")
            sys.exit()
    else:
        print("Wrong date format.")
        sys.exit()


if __name__ == "__main__":
    main()
