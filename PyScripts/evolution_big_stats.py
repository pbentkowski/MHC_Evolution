#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 02:38:00 2016

@author: piotr
"""
import re
# import sys
# import linecache as ln
import numpy as np
# import matplotlib.pyplot as plt
import bitstring as bts


def loadHostPopulation(FILE):
    '''Takes the file with the Host population HostGenomesFile.XXXX.csv and
    picks unique genes from it. Produces two lists: one containing ancestry of
    each gene (tags of all predecessors) and times when each mutation arose in
    the timeline.'''
    B_list = []
    Mut_tags = []
    Mut_times = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    continue
                else:
                    LL = line.split()
                    bb = bts.BitString(bin=LL[0]).int
                    if bb in B_list:
                        pass
                    else:
                        # print(LL[5::2])
                        B_list.append(bb)
                        tagz = LL[5::2]
                        tagz.append(LL[3])
                        Mut_tags.append(tagz)
                        timez = LL[4::2]
                        timez.append(LL[2])
                        Mut_times.append(timez)
        return Mut_tags, Mut_times
    except IOError as e:
        print("I/O error({0}) in".format(e.errno) +
              " loadTheHostPopulation(): {0}".format(e.strerror))


def loadRawBitstrings(FILE):
    """ """
    B_list = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    continue
                else:
                    LL = line.split()
                    bb = LL[0]
                    if bb in B_list:
                        pass
                    else:
                        B_list.append(bb)
        return B_list
    except IOError as e:
        print("I/O error({0}) in".format(e.errno) +
              " loadTheHostPopulation(): {0}".format(e.strerror))


def findTheOnesAtBeginning(Mut_tags, jj=0):
    """ """
    ll = []
    for itm in Mut_tags:
        try:
            if itm[jj] in ll:
                pass
            else:
                ll.append(itm[jj])
        except:
            pass
    return ll


def numberOfMutList(Mut_tags):
    """ """
    ll = []
    for itm in Mut_tags:
        ll.append(len(itm))
    return np.array(ll)


def findMRCA(Mut_tags, Mut_times):
    """Finds the tag, time stamp and index of the most recent common ancestor
    gene from the list of all genes at the population snapshot."""
    if len(findTheOnesAtBeginning(Mut_tags, 0)) != 1:
        print("The most recent common ancestor cannot be established.",
              "There is more than one ancestral gene at the root.")
        return None, np.nan, np.mnan
    mutNumb = numberOfMutList(Mut_tags)
    maxx = np.max(mutNumb)
    theMRCAtag = Mut_tags[0][0]
    ii = 0
    for x in range(maxx):
        if len(findTheOnesAtBeginning(Mut_tags, x)) == 1:
            theMRCAtag = findTheOnesAtBeginning(Mut_tags, x)[0]
            ii = x
        else:
            break
    return theMRCAtag, int(Mut_times[0][ii]), ii


def timeOfExistence(Mut_tags, Mut_times):
    """ """
    times = []
    for itm in Mut_times:
        for ii in range(1, len(itm)):
            times.append(int(itm[ii]) - int(itm[ii-1]))
        return np.array(times)
