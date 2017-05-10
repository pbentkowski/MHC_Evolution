#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:27:47 2017

@author: piotr
"""
import re
import linecache as ln
import numpy as np
import matplotlib.pyplot as plt


def getUniqueGenes(FILE, firstIndex=1):
    """ """
    uniqueGenes = []
    with open(FILE) as infile:
        for i, line in enumerate(infile):
            if re.search(r"#", line):
                continue
            elif(i < firstIndex):
                continue
            else:
                LL = line.split()
                for j in range(1, len(LL)):
                    geneID = int(LL[j])
                    if(geneID in uniqueGenes):
                        pass
                    else:
                        uniqueGenes.append(geneID)
    return uniqueGenes


def getIndexGivenGeneID(FILE, geneID, firstIndex=1):
    """ """
    indexInFile = []
    strID = str(geneID)
    with open(FILE) as infile:
        for i, line in enumerate(infile):
            if re.search(r"#", line):
                continue
            elif(i < firstIndex):
                continue
            else:
                if re.search(" " + strID, line):
                    LL = line.split()
                    indexInFile.append((i + 1, LL.index(strID)))
    return np.array(indexInFile, dtype=int)


def getDataOfGivenGeneByID(FILE, idxArr, firstIndex=1):
    """ """
    dataInFile = []
    for itm in idxArr:
        l = re.split(" ", ln.getline(FILE, itm[0]))
        dataInFile.append((int(l[itm[1]])))
    return dataInFile


#==============================================================================
# abundance_97 = np.array(ibs.getDataOfGivenGeneByID("InfectionGeneNumb.csv", idxArr_97), dtype=int)
# presented_97 = np.array(ibs.getDataOfGivenGeneByID("InfectionGeneNumbPatho.csv", idxArr_97), dtype=int)
#==============================================================================
