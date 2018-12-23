#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Your doc string here, please...


Created on Fri Feb  3 15:30:05 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import os
import sys
import json
import numpy as np
# import linecache as ln
import matplotlib.pyplot as plt
from scipy.stats import linregress
# depends on this packedge of mine:
import packed_plots_of_MHC_alleles as ppma


# """Data type for storing processed data"""
outType = np.dtype([('VAR', 'f8'), ('VARX', 'f8'), ('slope', 'f8'),
                    ('intercept', 'f8'), ('R2', 'f8'), ('p', 'f8'),
                    ('sd_err', 'f8'), ('pathoNumb', 'f8'),
                    ('sourceDir', 'S99')])
handyType = np.dtype([('VAR', 'f8'), ('VARX', 'f8'), ('slope', 'f8'),
                      ('intercept', 'f8'), ('R2', 'f8'), ('p', 'f8'),
                      ('sd_err', 'f8'), ('pathoNumb', 'f8')])


def importComputedData(dataFile):
    """ """
    try:
        dd = np.genfromtxt(dataFile, usecols=(0, 1, 2, 3, 4, 5, 6, 7),
                           skip_header=1, dtype=handyType)
        return dd
    except Exception:
        print("ERROR in importComputedData(): Cannot load the data from file:",
              dataFile)


def the_line(x, a, b):
    return a * x + b


def loadHostPopulation(FILE):
    '''Takes the file with all the hosts data loads it to a list. Each
    individual is loaded as a list of bit strings.And the population is a list
    of individuals.'''
    LL = []
    ll = []
    nextHost = False
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    if nextHost:
                        LL.append(ll)
                        ll = []
                else:
                    nextHost = True
                    ll.append(line.split()[0])
        LL.append(ll)
#        print j + 1
        return LL
    except IOError as e:
        print("\nI/O error({0}) in".format(e.errno),
              "loadTheHostPopulation(): {0}".format(e.strerror))
        return None


def loadPathoExposed(FILE):
    """ """
    LL = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    ll = line.split()
                    if int(ll[7]) > 0:
                        LL.append(ll[11:-1])
                    else:
                        LL.append([])
                else:
                    continue
        return LL
    except IOError as e:
        print("\nI/O error({0}) in".format(e.errno),
              "loadPathoExposed(): {0}".format(e.strerror))
        return None


def uniqueMhcInHostOnly(hostPopList):
    """ """
    uniqHosts = []
    for indv in hostPopList:
        ll = []
        for mhc in indv:
            if mhc not in ll:
                ll.append(mhc)
        uniqHosts.append(ll)
    return uniqHosts


def calculateTheNumbers(hostPopList, pathoExposed):
    """ """
    uniq = uniqueMhcInHostOnly(hostPopList)
    LL = []
    for indv in uniq:
        LL.append(len(indv))
    uniqNumb = np.array(LL)
    ll = []
    for ii in pathoExposed:
        ll.append(len(ii))
    pathoNumb = np.array(ll)
    return uniqNumb, pathoNumb


def plotMHCvsPathoPresent(uniqNumb, pathoNumb, slope, intercept,
                          dirr=os.getcwd(), jitter=0.05):
    """ """
    plt.cla()
    jitterX = jitter * np.random.randn(uniqNumb.shape[0])
    jitterY = jitter * np.random.randn(pathoNumb.shape[0])
#    slope, intercept, r_val, p_val, std_err = linregress(uniqNumb, pathoNumb)
#    print("R^2 =", r_val**2)
#    print("p =", p_val)
    FS = 14
    TicksFS = 12
    plt.figure(1, figsize=(16, 16))
    plt.scatter(uniqNumb + jitterX, pathoNumb + jitterY)
    plt.plot(uniqNumb, the_line(uniqNumb, slope, intercept), 'r-')
    plt.xlabel("number of unique MHC alleles in individual", fontsize=FS)
    plt.ylabel("number of infections presented\nto immune system", fontsize=FS)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
#    plt.grid(True)
#    plt.show()
    plt.savefig(dirr + "/infect_vs_MHCnumb.png")


def getTheData(theStartDate, templateList, dirr=os.getcwd()):
    """Walking the dir using Python 3.6. Variable theStartDate has to be
    a datetime.date() data type."""
    vv = ppma.lookForVARinList(templateList)
    datOut = []
    dataOrdering = ['VAR', 'VARX', 'slope', 'intercept']
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.json') and
               ppma.loadTheDateFromParamFile(filepath) >= theStartDate):
                paramzList = ppma.loadParamSettings(filepath)
                if ppma.compareParams(templateList, paramzList):
                    with open(filepath) as f:
                        prms = json.load(f)
                    path_spp = float(prms['number_of_pathogen_species'])
                    lg = prms['number_of_host_generations']
                    genomeFileName = "HostGenomesFile." + str(lg) + ".csv"
                    genomeFileName = os.path.join(dirName, genomeFileName)
#                    print(genomeFileName)
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    try:
                        print(dirName, end=' : ')
                        pathos = loadPathoExposed(genomeFileName)
                        hosts = loadHostPopulation(genomeFileName)
                        if hosts is None or pathos is None:
                            print("Failed to read data")
                            continue
                        else:
                            print("Done")
                    except Exception:
                        print("ERROR in getTheData(): cant's load the host",
                              "population data")
                        continue
                    uniqNumb, pathoNumb = calculateTheNumbers(hosts, pathos)
                    uniqNumb = np.hstack((uniqNumb, 0))
                    pathoNumb = np.hstack((pathoNumb, 0))
                    # slope, intercept, r_val, p_val, std_err
                    data = linregress(uniqNumb, pathoNumb)
                    plotMHCvsPathoPresent(uniqNumb, pathoNumb,
                                          data[0], data[1], dirName)
                    datOut.append((var, varx, data[0], data[1], data[2]**2,
                                   data[3], data[4], path_spp, dirName))
    datOut = np.array(datOut, dtype=outType)
    return np.sort(datOut, order=dataOrdering)


def plotBoxesSlopes(handyArr, xlabels="The values", yMax=0):
    """ """
    if handyArr.dtype == handyType:
        spp = np.unique(handyArr['VARX'])
    else:
        print("ERROR in plotHistograms(): wrong numpy data type. It should",
              "be:", handyType)
        return None
    lbls = spp.astype(int)
    ll = []
    fs = 16
    tkfs = 14
    plt.figure(figsize=(9, 6))
    for itm in spp:
        ll.append(handyArr[handyArr['VARX'] == itm]['slope'])
    boxprops = dict(linestyle='-', linewidth=2.5, color='k')
    medianprops = dict(linestyle='-', linewidth=2.5)
    whiskerprops = dict(linewidth=2.5)
    capprops = dict(linewidth=2.5)
    flierprops = dict(markersize=10)
    plt.boxplot(ll, labels=lbls, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops)
    plt.xlabel(xlabels, fontsize=fs)
    plt.ylabel("linear regression slope value", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    if yMax > 0:
        plt.ylim(ymax=yMax)
    plt.grid(axis='y')
    plt.show()


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 3:
        print("Two arguments are needed:")
        print("  1. Give a starting date. It has to be in yyyy-mm-dd format.")
        print("  2. Give the path to template file.")
        print("  3. Give the name of the output file.")
        sys.exit()
    startDate = None
    headerr = 'VAR VARX slope intercept R2 p_value sr_err patho_number '\
        + 'sourceDir'
    try:
        startDate = ppma.readDate(sys.argv[1])
        outputFile = str(sys.argv[3])
    except ValueError:
        print("Cannot convert argument #1 to a date format.")
        sys.exit()
    if startDate:
        try:
            template = ppma.loadParamSettings(sys.argv[2])
#            x_Label = ppma.getVarxLabel(sys.argv[2])
        except Exception:
            print("Cannot load the template file. Exiting.")
            sys.exit()
        try:
            print("Computing data...")
            theData = getTheData(startDate, template)
        except Exception:
            print("Failed to process the data. Some serious issues arose.")
            sys.exit()
        if len(theData):
            FMT = '%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %s'
            open(outputFile, 'w').close()
            np.savetxt(outputFile, theData, fmt=FMT, header=headerr,
                       comments='#')
            for itm in theData:
                for ii in range(len(itm) - 1):
                    print(itm[ii], "\t", end=" ")
                print()
            print("Check the output file:", str(os.getcwd()) + "/" +
                  outputFile + " for details.")
        else:
            print("No data files matching the criterions were found.",
                  "Specify your template file.")
            sys.exit()
    else:
        print("Wrong date format.")
        sys.exit()


if __name__ == "__main__":
    main()
