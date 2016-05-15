#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Walks the directory tree looking for model runs which are characterized by same
parametrisation as the template file provided by the user. Then process these
results by fancy stats and plot that processed output on a nice plot.

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


# """Data type for loading data from files HostsGeneDivers.csv"""
datType = np.dtype([('time', np.int), ('pop_size', np.int),
                    ('tot_num_of_genes', np.int), ('num_of_MHC_types', np.int),
                    ('Shannon_indx', np.float), ('mean_fitness', np.float),
                    ('std_fitness', np.float)])
outType = np.dtype([('VAR', 'f8'), ('VARX', 'f8'), ('meanAllel', 'f8'),
                    ('stdAllel', 'f8'), ('slope', 'f8')])


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
        for ii, itm in enumerate(zip(template[::-3], paramz[::-3])):
            try:
                ITM_0 = float(itm[0])
                ITM_1 = float(itm[1])
            except:
                ITM_0 = str(itm[0])
                ITM_1 = str(itm[1])
            if itm[0] == "VARX" or itm[0] == "VAR" or ii <= 1:
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
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.csv') and
               loadTheDateFromParamFile(filepath) >= theStartDate):
                paramzList = loadParamSettings(filepath)
                if compareParams(template, paramzList):
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    dataFilePath = os.path.join(dirName, "HostsGeneDivers.csv")
                    data = np.genfromtxt(dataFilePath, dtype=datType)
                    c0, c1 = poly.polyfit(data['time'][EqPt::],
                                          data['num_of_MHC_types'][EqPt::], 1)
                    meanAlle = data['num_of_MHC_types'][EqPt::].mean()
                    stdAlle = data['num_of_MHC_types'][EqPt::].std()
                    datOut.append((var, varx, meanAlle, stdAlle, c1))
    datOut = np.array(datOut, dtype=outType)
    return np.sort(datOut,
                   order=['VAR', 'VARX','meanAllel', 'stdAllel', 'slope'])


def main():
    """ """
    if len(sys.argv) <= 2:
        print("Two arguments are needed:")
        print("  1. Give a starting date. It has to be in yyyy-mm-dd format.")
        print("  2. Give the path to template file.")
        sys.exit()
    try:
        startDate = readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert argument #1 to a date format.")
        sys.exit()
    if startDate:
        try:
            template = loadParamSettings(sys.argv[2])
        except:
            print("Cannot load the template file. Exiting.")
            sys.exit()
        theData = getTheData(startDate, template)
        for itm in theData:
            for ii in itm:
                print(ii, "\t", end=" ")
            print()
    else:
        print("Wrong date format.")
        sys.exit()


if __name__ == "__main__":
    main()
