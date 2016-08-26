#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Calculates if MHC genes are correlated in co-occurrence in individuals as
a result e.g. sex selection

Created on Fri Aug 19 18:22:31 2016
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
# import sys
import re
# import linecache as ln
import numpy as np
import matplotlib.pyplot as plt
import bitstring as bts


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
                    bitt = line.split()[0]
                    ll.append(bts.BitString(bin=bitt).int)
        LL.append(np.array(ll, dtype=int))
#        print j + 1
        return LL
    except IOError as e:
        print("I/O error({0}) in".format(e.errno),
              "loadHostPopulation(): {0}".format(e.strerror))


def calculateCooccurMatrix(hostList):
    """Takes the list of hosts and their MHC genes and calculates co-occurrence
    of genes found in the population."""
    if len(hostList) == 0:
        print("ERROR in calculateCooccurMatrix(): the hosts list is empty.")
        return None
    # === fetch unique MHC genes ===
    geneList = []
    for indv in hostList:
        for gene in indv:
            if gene in geneList:
                pass
            else:
                geneList.append(gene)
    geneVec = np.array(geneList, dtype=int)
    N = geneVec.shape[0]
    # === fetch frequencies of occurrence of MHC genes ===
    geneFreq = np.zeros(N)
    for k, gg in enumerate(geneVec):
        cc = 0
        for indv in hostList:
            if(gg in indv):
                cc += 1
        geneFreq[k] = cc
    # === calculate MHC genes co-occurrence matrix ===
    geneFreqMtx = np.zeros((N, N))
    cooccMtx = np.zeros((N, N))
    for i, one in enumerate(geneVec):
        for j, two in enumerate(geneVec):
            count = 0
            for indv in hostList:
                if(one in indv and two in indv):
                    count += 1
            cooccMtx[i, j] = count
            geneFreqMtx[i, j] = np.max([geneFreq[i], geneFreq[j]])
    # === return normalized co-occurrence matrix ===
    return geneFreq, cooccMtx / geneFreqMtx
