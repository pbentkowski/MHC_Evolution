#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
                The doc string here, please....

Created on Fri Nov  3 19:19:04 2017

@author: Piotr Bentkowski :: bentkowski.piotr@gmail.com
"""
import os
import sys
import re
import linecache as ln
import infectionStatsBig as isb


def getTheData(firstGen, numbOfRand, minLastTime, avgWay='median',
               dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type."""
    datOut = []
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.csv')):
                print("processing data in:", dirName)
                geneList = isb.pickRandomMHCs(dirName + "/InfectionGeneID.csv",
                                              firstGen, numbOfRand)
                paramsFile = os.path.join(dirName, 'InputParameters.csv')
                ll = re.split(" ", ln.getline(paramsFile, 7))
                hostPopSize = float(ll[2])
                ll = re.split(" ", ln.getline(paramsFile, 9))
                pathoSppNum = float(ll[2])
                ll = re.split(" ", ln.getline(paramsFile, 20))
                alphaFact = float(ll[2])
                mhcStats = isb.calcRelatFittManyMHC(geneList, hostPopSize,
                                                    avgWay, dirName)
                mhcStats, procFixed = isb.removeShortLivedMHC(mhcStats,
                                                              minLastTime)
                avgImmCompts = isb.calcAverageForOneRun(mhcStats, avgWay)
                datOut.append((alphaFact, pathoSppNum, procFixed,
                               avgImmCompts))
    return datOut


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 3:
        print("Three arguments are needed:")
        print("  1. Host generation after which stats will be obtained.")
        print("  2. Number of genes that need to be selected for analysis.")
        print("  3. Minimum number of generations an MHC must exist.")
        sys.exit()
#    startDate = None
#    headerr = 'VAR VARX meanAllel stdAllel slope indvMean indvSTD meanFitt '\
#        + 'stdFitt meanCvFitt stdCvFitt sourceDir'
    if 1:
        # third argument is very important
        theData = getTheData(int(sys.argv[1]), int(sys.argv[2]),
                             int(sys.argv[3]))
        print(theData)
#    except Exception:
#        print("Failed to process the data. Some serious issues arose.",
#              "Check if the cut-off host generation for calculating stats",
#              "is smaller than the total number of host generations.")
#        sys.exit()


if __name__ == "__main__":
    main()
