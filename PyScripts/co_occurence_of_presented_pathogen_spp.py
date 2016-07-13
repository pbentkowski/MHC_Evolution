#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Your doc string goes here...

Created on Wed Jul 13 18:26:33 2016
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import os
import linecache as ln
import numpy as np


def loadParams(dirname):
    """ """
    try:
        paramsFile = os.path.join(dirname, 'InputParameters.csv')
        l = re.split(" ", ln.getline(paramsFile, 8))
        popSize = int(l[2])
        l = re.split(" ", ln.getline(paramsFile, 9))
        sppNumber = int(l[2])
        return popSize, sppNumber
    except:
        print("ERROR in loadParams(): cannot find parameters in file")
        return None, None


def loadPresentedSpecies(filepath):
    """ """
    LL = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if re.search(" ===", line):
                    LL.append(list(map(int, line.split("are:")[1].
                              split(" ===")[0].split())))
                else:
                    pass
        return LL
    except:
        print("ERROR in loadPresentedSpecies(): Cannot load the presented",
              "pathogen species.")
        return None


def calculateCoocurenceMtx(LL, sppNumber, popSize):
    """ """
    mtx = np.zeros((sppNumber, sppNumber))
    for i in range(sppNumber):
        for j in range(i+1, sppNumber):
            for itm in LL:
                if i in itm and j in itm:
                    mtx[i, j] += 1
                    mtx[j, i] += 1
    return mtx / float(popSize)
