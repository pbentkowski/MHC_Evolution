#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Iterates through the directories looking for the data regarding the infection
data statistics. Generates a file containing a list of stats from each
individual run (one per line) with its main parameters.

Created on Fri Nov  3 19:19:04 2017

@author: Piotr Bentkowski :: bentkowski.piotr@gmail.com
"""
import os
import sys
import re
import linecache as ln
import infectionStatsBig as isb


def getTheData(firstGen, numbOfRand, minLastTime, avgWay='mean',
               dirr=os.getcwd()):
    """Walking the dir in search of data using the os.walk() mechanism"""
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
                mhcStats, totGen, genFix = isb.removeShortLivedMHC(mhcStats,
                                                                   minLastTime)
                avgImmCompts = isb.calcAverageForOneRun(mhcStats, avgWay)
                datOut.append((alphaFact, pathoSppNum, totGen, genFix,
                               avgImmCompts))
                print
    return datOut


def argInfo():
    print("Three arguments are needed:")
    print("  1. Host generation after which stats will be obtained.")
    print("  2. Number of genes that need to be selected for analysis.")
    print("  3. Minimum number of generations an MHC must exist.")
    print("  4. How to average the data? 'mean' or 'median'?")


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 4:
        argInfo()
        sys.exit()
    file_c = "#alpha_factor numb_of_patho_spp genesTested genesFixed " \
        + "array_with_avg_immunocompetence\n"
    if sys.argv[4] == 'mean' or sys.argv[4] == 'median':
        pass
    else:
        print("Averaging methods must be 'mean' or 'median'!")
        sys.exit()
    try:
        generStart = int(sys.argv[1])
        numbGenes = int(sys.argv[2])
        minGeneAge = int(sys.argv[3])
        f_name = "immuno_" + sys.argv[1] + "_" + sys.argv[2] + "_" \
                 + sys.argv[3] + "_" + sys.argv[4] + ".dat"
    except Exception:
        print("Cannot convert arguments to integers. Check your params")
        argInfo()
        sys.exit()
    try:
        theData = getTheData(generStart, numbGenes, minGeneAge, 'median')
    except Exception:
        print("Failed to process the data. Some serious issues arose.",
              "Check if the cut-off host generation for calculating stats",
              "is smaller than the total number of host generations.")
        sys.exit()
    for itm in theData:
        file_c += str(itm[0]) + ',' + str(itm[1]) + ',' + str(itm[2]) + ',' \
            + str(itm[3]) + ','
        for ii in itm[4]:
            file_c += "%6.3f;" % ii
        file_c += "\n"
    dataFile = open(f_name, 'w')
    dataFile.write(file_c)
    dataFile.close
    print("DONE! Check file", f_name, "for data output.")


if __name__ == "__main__":
    main()
