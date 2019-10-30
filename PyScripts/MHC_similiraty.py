#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Calculates and plots similarities between antigens and genes in pathogen and
host populations. Does it separately for hosts and for pathogens.

Created on Tue Nov  3 16:10:58 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import re
import json
# import linecache as ln
import numpy as np
import matplotlib.pyplot as plt


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
    except Exception:
        print("ERROR in anti_gen_similiraty.bitSimWhinIndiv():",
              "Can't load the data!")
        return np.NaN


def hamDistWhinIndiv(BitLyst):
    """Hamming distances between genes in one pathogen."""
    try:
        N = len(BitLyst)
        compArr = np.ones(int((N*(N-1))/2))
        k = 0
        for ii, itmOne in enumerate(BitLyst):
            for jj in np.arange(ii+1, N, 1):
                compArr[k] = hamming_distance(itmOne, BitLyst[jj])
                k += 1
        return compArr
    except Exception:
        print("ERROR in anti_gen_similiraty.bitSimWhinIndiv():",
              "Can't load the data!")
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
    except Exception:
        print("ERROR in anti_gen_similiraty.bitSimBetweenIndv():",
              "Can't proccess the data!")
        return np.NaN


def hamDistBetweenIndv(indOne, indTwo):
    """Takes two sets of antigens (tow individual pathogens) and compares them
    (antigen by antigen) according to Hamming distance."""
    try:
        compArr = np.zeros(len(indOne) * len(indTwo))
        k = 0
        for itmOne in indOne:
            for itmTwo in indTwo:
                compArr[k] = hamming_distance(itmOne, itmTwo)
                k += 1
        return compArr
    except Exception:
        print("ERROR in anti_gen_similiraty.hamDistBetweenIndv():",
              "Can't proccess the data!")
        return np.NaN


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
        print("I/O error({0}) in".format(e.errno),
              "loadTheHostPopulation(): {0}".format(e.strerror))


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
        print("I/O error({0}) in".format(e.errno),
              "loadPathoExposed(): {0}".format(e.strerror))


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


def plotMHCvsPathoPresent(hostPopList, pathoExposed):
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
    plt.figure(1, figsize=(12, 8))
    plt.plot(uniqNumb, pathoNumb, 'o')
    plt.grid(True)
    plt.show()


def bitSimAll(Popul, simmes=7):
    """Calculates bit string similarity for all individuals (mean is calculated
    for a single Pathogen. """
    DD = []
    for indv in Popul:
        DD.append(np.mean(bitSimWhinIndiv(indv, simmes)))
    return np.array(DD)


def bitSimInterIndv(Popul, simms):
    """ """
    N = len(Popul)
    compArr = np.zeros(2*N)
    rndIndx_1 = np.random.randint(0, N, 2*N)
    rndIndx_2 = np.random.randint(0, N, 2*N)
    for kk, comp in enumerate(compArr):
        while (rndIndx_1[kk] == rndIndx_2[kk]):
            rndIndx_1[kk] = np.random.randint(0, N)
        compArr[kk] = np.mean(bitSimBetweenIndv(Popul[rndIndx_1[kk]],
                              Popul[rndIndx_2[kk]], simms))
    return compArr


def hamDisthistAll(Popul):
    """Mean Hamming distances for all individuals (mean is calculated for
    a single Pathogen."""
    DD = []
    for indv in Popul:
        DD.append(np.mean(hamDistWhinIndiv(indv)))
    return np.array(DD)


def hamDistInterIndv(Popul):
    """ """
    N = len(Popul)
    compArr = np.zeros(2*N)
    rndIndx_1 = np.random.randint(0, N, 2*N)
    rndIndx_2 = np.random.randint(0, N, 2*N)
    for kk, comp in enumerate(compArr):
        while (rndIndx_1[kk] == rndIndx_2[kk]):
            rndIndx_1[kk] = np.random.randint(0, N)
        compArr[kk] = np.mean(hamDistBetweenIndv(Popul[rndIndx_1[kk]],
                              Popul[rndIndx_2[kk]]))
    return compArr


def main():
    if len(sys.argv) <= 2:
        print("Give the names of two files with data. One at the begging of",
              "simulation and the second at the end.")
        sys.exit()
    try:
        with open('InputParameters.json') as f:
            prms = json.load(f)
        pathSpp = prms['number_of_pathogen_species']
        bitfit = prms['number_of_bits_per_antigen']
        print("No. of pathogen species =", int(pathSpp))
        MHC_len = int(prms['number_of_bits_per_gene'])
        print("Length of the bit string =", MHC_len)
        print("Length of bit string fit =", bitfit)
    except Exception:
        print("Can't find parameter file! You may be in a wrong directory.")
        sys.exit()
    try:
        L_init = loadHostPopulation(sys.argv[1])
        print("First file loaded!")
    except Exception:
        print("Can't load file named", sys.argv[1], ". Check if it exists.")
        sys.exit()
    try:
        L_endd = loadHostPopulation(sys.argv[2])
        print("Second file loaded!")
    except Exception:
        print("Can't load file named", sys.argv[2], ". Check if it exists.")
        sys.exit()
#    F_init = bitSimAll(L_init, bitfit)
    F_init = hamDisthistAll(L_init)
    F_init = F_init[~np.isnan(F_init)]
    print("Within genome similarities in the First file have been calculated!")
#    F_endd = bitSimAll(L_endd, bitfit)
    F_endd = hamDisthistAll(L_endd)
    F_endd = F_endd[~np.isnan(F_endd)]
    print("Within genome similarities in the Second file have been" +
          " calculated!")
#    E_init = bitSimInterIndv(L_init, bitfit)
    E_init = hamDistInterIndv(L_init)
    E_init = E_init[~np.isnan(E_init)]
    print("Between individual similarities in the First file have",
          "been calculated!")
#    E_endd = bitSimInterIndv(L_endd, bitfit)
    E_endd = hamDistInterIndv(L_endd)
    E_endd = E_endd[~np.isnan(E_endd)]
    print("Between individual similarities in the Second file have",
          "been calculated!")
    # === More generic plot ===
    ax_label = 20
    T_label = 24
    TicksFS = 18
    transs = 0.8
    plt.figure(1, figsize=(16, 8))
    plt.subplot(121)
    plt.hist(F_init, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("Start of simulation", fontsize=T_label)
    plt.xlabel("Within genome similarity measure", fontsize=ax_label)
    plt.ylabel("Frequency of occurrence", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
    plt.xlim(0, MHC_len)
    plt.subplot(122)
    plt.hist(F_endd, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("End of simulation", fontsize=T_label)
    plt.xlabel("Within genome similarity measure", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
    plt.xlim(0, MHC_len)
    plt.savefig("HOST_sim_one.png")
    #  === Now the detailed plot! ===
    plt.figure(2, figsize=(16, 8))
    plt.subplot(121)
    plt.hist(E_init, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("Start of simulation", fontsize=T_label)
    plt.xlabel("Inter-individual similarity measure", fontsize=ax_label)
    plt.ylabel("Frequency of occurrence", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
    plt.xlim(0, MHC_len)
#    plt.ylim(ymax=200)
    plt.subplot(122)
    plt.hist(E_endd, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("End of simulation", fontsize=T_label)
    plt.xlabel("Inter-individual similarity measure", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
    plt.xlim(0, MHC_len)
#    plt.ylim(ymax=200)
    plt.savefig("HOST_sim_two.png")

    plt.show()

    print("DONE!")


if __name__ == "__main__":
    main()
