#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Fri Oct 23 17:05:28 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import re
import linecache as ln
import numpy as np
#from bitstring import BitArray


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def bitSimInRow(s1, s2, sim_measure):
    """Compares two bit string of the same length. They need to have a defined
    number of common bits IN A ROW to be called \"similar\"."""
    if(len(s1) == len(s2)):
        couter = 0
        for itm in zip(s1, s2):
            if(itm[0] == itm[1]):
                couter += 1
            else:
                couter = 0
            if(couter >= sim_measure):
                return True
        if(couter < sim_measure):
            return False
    else:
        return False


def bitSimWhinIndiv(BitLyst, sim_measure=7):
    """ """
    try:
        N = len(BitLyst)
        compArr = np.ones((N*(N-1))/2)
        k = 0
        for ii, itmOne in enumerate(BitLyst):
            for jj in np.arange(ii+1, N, 1):
                if bitSimInRow(itmOne, BitLyst[jj], sim_measure):
                    compArr[k] = 1
                else:
                    compArr[k] = 0
                k += 1
        return np.mean(compArr)
    except:
        print "ERROR in anti_gen_similiraty.bitSimWhinIndiv():",
        print "Can't load the data!"
        return np.NaN


def hamDistWhinIndiv(BitLyst):
    """Hamming distances between genes in one pathogen."""
    try:
        N = len(BitLyst)
        compArr = np.ones((N*(N-1))/2)
        k = 0
        for ii, itmOne in enumerate(BitLyst):
            for jj in np.arange(ii+1, N, 1):
                compArr[k] = hamming_distance(itmOne, BitLyst[jj])
                k += 1
        return compArr
    except:
        print "ERROR in anti_gen_similiraty.bitSimWhinIndiv():",
        print "Can't load the data!"
        return np.NaN


def bitSimBetweenIndv(indOne, indTwo, sim_measure=7):
    """ """
    try:
        compArr = np.zeros(len(indOne))
        for k, itmOne in enumerate(indOne):
            for itmTwo in indTwo:
                if bitSimInRow(itmOne, itmTwo, sim_measure):
                    compArr[k] = 1.
        return compArr
    except:
        print "ERROR in anti_gen_similiraty.bitSimBetweenIndv():",
        print "Can't proccess the data!"
        return np.NaN


def loadThePopulation(FILE):
    '''Takes the file with all the pathogen data loads it to a list dividing
    the population into species and individuals.'''
    LL = []
    spp_list = []
    nextPatho = False
    ll = []
    l = re.split(" ", ln.getline(FILE, 3))
    old_spec = int(l[5])
    try:
        with open(FILE) as infile:
            for j, line in enumerate(infile):
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    try:
                        new_spec = int(line.split()[4])
                    except:
                        pass
                    if nextPatho:
                        if new_spec == old_spec:
                            spp_list.append(ll)
                            ll = []
                            endOfSpp = False
                        else:
                            if endOfSpp:
                                LL.append(spp_list)
                                spp_list = []
                                spp_list.append(ll)
                                ll = []
                                old_spec = new_spec
#                                print j + 1
                            else:
                                spp_list.append(ll)
                                ll = []
                            endOfSpp = True
                    nextPatho = True
                else:
                    ll.append(line.split()[0])
        spp_list.append(ll)
        LL.append(spp_list)
#        print j + 1
        return LL
    except IOError as e:
        print "I/O error({0}) in loadThePopulation(): {1}".format(e.errno,
                                                                  e.strerror)


def hamDisthistAll(LL):
    """Mean Hamming distances for all individuals (mean is calculated for
    a single Pathogen."""
    DD = []
    for sp in LL:
        for indv in sp:
            DD.append(np.mean(hamDistWhinIndiv(indv)))
    return np.array(DD)


def main():
    if len(sys.argv) <= 1:
        print "Give the name of the file with data."
        sys.exit()
    loadThePopulation(sys.argv[1])
    print "DONE!"


if __name__ == "__main__":
    main()
