#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Your doc string goes here...

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


def loadParamSettings(filepath):
    """ """
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
    """ """
    same = False
    for ii, itm in enumerate(zip(template, paramz)):
        try:
            ITM_0 = float(itm[0])
            ITM_1 = float(itm[1])
        except:
            ITM_0 = str(itm[0])
            ITM_1 = str(itm[1])
        if itm[0] == "VARX" or itm[0] == "VAR" or ii == 0:
            print(itm)
        elif ITM_0 == ITM_1:
            same = True
        else:
            same = False
            break
    return same


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


def printTheParams3(theStartDate, dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type."""
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.csv') and
               loadTheDateFromParamFile(filepath) >= theStartDate):
                strr = ""
                with open(filepath, 'r') as f:
                    for ii, line in enumerate(f):
                        l = line.split()
                        if re.search("#", line):
                            if re.search("# Other_information:", line):
                                break
                            else:
                                continue
                        elif ii > 2 and len(l) > 1:
                            if l[2] == "YES":
                                strr += "10 "
                            elif l[2] == "NO":
                                strr += "11 "
                            else:
                                strr += l[2] + " "
                print(strr)


def main():
    """ """
    if len(sys.argv) <= 1:
        print("Give a starting date. It has to be in yyyy-mm-dd format.")
        sys.exit()
    try:
        startDate = readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert the input do a date format.")
        sys.exit()
    if startDate:
        printTheParams3(startDate)
    else:
        print("Wrong date format.")
        sys.exit()

if __name__ == "__main__":
    main()
