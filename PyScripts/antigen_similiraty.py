#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Checks the similarities in all the antigens in the hole pathogen population.
Function loadThePopulation(FILE) loads the file created by the modelling
framework into a handy data structure and the rest calculates various
statistics.

Created on Fri Oct 23 17:05:28 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import re
import linecache as ln
import numpy as np
import matplotlib.pyplot as plt
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
    """Calculates similarity between antigens within individual by comparing
    each antigen pairwise with all the others."""
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
    """Takes two sets of antigens (two individual pathogens) and compares them
    (antigen by antigen) according to their fit to MHC"""
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


def hamDistBetweenIndv(indOne, indTwo):
    """Takes two sets of antigens (two individual pathogens) and compares them
    (antigen by antigen) according to Hamming distance."""
    try:
        compArr = np.zeros(len(indOne) * len(indTwo))
        k = 0
        for itmOne in indOne:
            for itmTwo in indTwo:
                compArr[k] = hamming_distance(itmOne, itmTwo)
                k += 1
        return compArr
    except:
        print "ERROR in anti_gen_similiraty.hamDistBetweenIndv():",
        print "Can't proccess the data!"
        return np.NaN


def loadThePopulation(FILE):
    '''Takes the file with all the pathogen data loads it to a list dividing
    the population into species and individuals. Each individual is loaded as
    a list of bit strings. Each species is a list of individuals. And the
    population is a list of species lists.'''
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
        print "I/O error({0}) in".format(e.errno),
        print "loadTheHostPopulation(): {0}".format(e.strerror)


def bitSimAll(Popul, simmes=7):
    """Calculates bit string similarity for all individuals (mean is calculated
    for a single Pathogen. """
    DD = []
    for sp in Popul:
        for indv in sp:
            DD.append(np.mean(bitSimWhinIndiv(indv, simmes)))
    return np.array(DD)


def bitSimInsideSpec(Spec, simm):
    """Calculates bit string similarity between individuals within the species
    for a single species. Returns an array of individual-to-individual
    comparisons."""
    N = len(Spec)
    compArr = np.zeros(2*len(Spec))
    rndIndx_1 = np.random.randint(0, N, 2*N)
    rndIndx_2 = np.random.randint(0, N, 2*N)
    for kk, comp in enumerate(compArr):
        while (rndIndx_1[kk] == rndIndx_2[kk]):
            rndIndx_1[kk] = np.random.randint(0, N)
        compArr[kk] = np.mean(bitSimBetweenIndv(Spec[rndIndx_1[kk]],
                              Spec[rndIndx_2[kk]], simm))
    return compArr


def bitSimInterSpec(Popul, simm):
    """Calculates bit string similarity between individuals within the species
    for each species in the population. Return an NumPy array of mean
    similarities for species individually."""
    CMP = []
    for ii, spp1 in enumerate(Popul):
        for jj in np.arange(ii+1, len(Popul)):
            compArr = np.zeros(2*len(spp1))
            rndIndx_1 = np.random.randint(0, len(spp1), 2*len(spp1))
            rndIndx_2 = np.random.randint(0, len(Popul[jj]), 2*len(spp1))
            for kk, comp in enumerate(compArr):
                compArr[kk] = np.mean(bitSimBetweenIndv(spp1[rndIndx_1[kk]],
                                      Popul[jj][rndIndx_2[kk]]))
            CMP.append(compArr)
    return CMP


def hamDisthistAll(Popul):
    """Mean Hamming distances for all individuals (mean is calculated for
    a single Pathogen."""
    DD = []
    for sp in Popul:
        for indv in sp:
            DD.append(np.mean(hamDistWhinIndiv(indv)))
    return np.array(DD)


def hamDistInsideSpec(Spec):
    """Calculates mean Hamming distances for individuals within the species.
    Returns a NumPy array of mean Hamming distances for number of pairwise
    comparisons between individuals."""
    N = len(Spec)
    compArr = np.zeros(2*len(Spec))
    rndIndx_1 = np.random.randint(0, N, 2*N)
    rndIndx_2 = np.random.randint(0, N, 2*N)
    for kk, comp in enumerate(compArr):
        while (rndIndx_1[kk] == rndIndx_2[kk]):
            rndIndx_1[kk] = np.random.randint(0, N)
        compArr[kk] = np.mean(hamDistBetweenIndv(Spec[rndIndx_1[kk]],
                              Spec[rndIndx_2[kk]]))
    return compArr


def hamDistInterSpecies(Popul):
    """Calculates mean Hamming distances for individuals between different
    species. Returns a list of NumPy arrays with comparison of each two species
    (makes a number of comparisons between a pair of randomly selected
    individuals from two species)."""
    CMP = []
    for ii, spp1 in enumerate(Popul):
        for jj in np.arange(ii+1, len(Popul)):
            compArr = np.zeros(2*len(spp1))
            rndIndx_1 = np.random.randint(0, len(spp1), 2*len(spp1))
            rndIndx_2 = np.random.randint(0, len(Popul[jj]), 2*len(spp1))
            for kk, comp in enumerate(compArr):
                compArr[kk] = np.mean(hamDistBetweenIndv(spp1[rndIndx_1[kk]],
                                      Popul[jj][rndIndx_2[kk]]))
            CMP.append(compArr)
    return CMP


def main():
    if len(sys.argv) <= 2:
        print "Give the names of two files with data. One at the begging of",
        print "simulation and the second at the end."
        sys.exit()
    try:
        l = re.split(" ", ln.getline("InputParameters.csv", 9))
        spp_num = int(l[2])
        l2 = re.split(" ", ln.getline("InputParameters.csv", 6))
        bitfit = int(l2[2])
        print "No. of pathogen species =", spp_num
        print "Length of the bit string =",
        print int(re.split(" ", ln.getline("InputParameters.csv", 5))[2])
        print "Length of bit string fit =", bitfit
    except:
        print "Can't load the parameter file! You may be in a wrong directory."
        sys.exit()
    try:
        L_init = loadThePopulation(sys.argv[1])
        print "First file loaded!"
    except:
        print "Can't load file named", sys.argv[1], ". Check if it exists."
        sys.exit()
    try:
        L_endd = loadThePopulation(sys.argv[2])
        print "Second file loaded!"
    except:
        print "Can't load file named", sys.argv[2], ". Check if it exists."
        sys.exit()
    F_init = hamDistInterSpecies(L_init)
    print "Similarities in the First file have been calculated!"
    F_endd = hamDistInterSpecies(L_endd)
    print "Similarities in the Second file have been calculated!"
    # === More generic plot ===
    ax_label = 20
    T_label = 24
    TicksFS = 18
    transs = 0.8
    plt.figure(1, figsize=(16, 8))
    plt.subplot(121)
    xx = np.zeros(len(F_init))
    for ii, itm in enumerate(F_init):
        xx[ii] = np.mean(itm)
    plt.hist(xx, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("Start of simulation", fontsize=T_label)
    plt.xlabel("Inter-species similarity measure", fontsize=ax_label)
    plt.ylabel("Frequency of occurrence", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
#    plt.xlim(0., 1.)
    plt.subplot(122)
    xx = np.zeros(len(F_endd))
    for ii, itm in enumerate(F_endd):
        xx[ii] = np.mean(itm)
    plt.hist(xx, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("End of simulation", fontsize=T_label)
    plt.xlabel("Inter-species similarity measure", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
#    plt.xlim(0., 1.)
    plt.savefig("SPP_sim_one.png")

    #  === Now the detailed plot! ===
    plt.figure(2, figsize=(16, 8))
    plt.subplot(121)
    plt.hist(F_init, edgecolor="none")
    plt.title("Start of simulation", fontsize=T_label)
    plt.xlabel("Inter-species similarity measure", fontsize=ax_label)
    plt.ylabel("Frequency of occurrence", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
#    plt.xlim(0., 1.)
#    plt.ylim(ymax=200)
    plt.subplot(122)
    plt.hist(F_endd, edgecolor="none")
    plt.title("End of simulation", fontsize=T_label)
    plt.xlabel("Inter-species similarity measure", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
#    plt.xlim(0., 1.)
#    plt.ylim(ymax=200)
    plt.savefig("SPP_sim_two.png")

    #  === Now the detailed plot! ===
    divv = int(np.ceil(np.sqrt(spp_num)))
    plt.figure(3, figsize=(24, 20))
#    ax_label_2 = 10
    TicksFS_2 = 11
    for ii in xrange(spp_num):
        plt.subplot(divv, divv, ii)
        plt.hist(hamDistInsideSpec(L_endd[ii]), color=(0.3, 0.3, 0.3, transs),
                 edgecolor="none")
#        plt.xlabel("Within-species similarity measure", fontsize=ax_label_2)
#        plt.ylabel("Frequency of occurrence", fontsize=ax_label_2)
        plt.xticks(fontsize=TicksFS_2)
        plt.yticks(fontsize=TicksFS_2)
        plt.grid(True)
    plt.savefig("SPP_sin_within_stop.png")

    plt.figure(4, figsize=(24, 20))
#    ax_label_2 = 10
    for ii in xrange(spp_num):
        plt.subplot(divv, divv, ii)
        plt.hist(hamDistInsideSpec(L_init[ii]), color=(0.3, 0.3, 0.3, transs),
                 edgecolor="none")
#        plt.xlabel("Within-species similarity measure", fontsize=ax_label_2)
#        plt.ylabel("Frequency of occurrence", fontsize=ax_label_2)
        plt.xticks(fontsize=TicksFS_2)
        plt.yticks(fontsize=TicksFS_2)
        plt.grid(True)
    plt.savefig("SPP_sin_within_init.png")

    plt.show()

    print "DONE!"


if __name__ == "__main__":
    main()
