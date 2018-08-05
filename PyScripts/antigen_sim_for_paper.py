#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 00:16:29 2018

@author: piotr
"""
import re
import sys
import linecache as ln
import numpy as np
import matplotlib.pyplot as plt
import antigen_similiraty as asim


def main():
    if len(sys.argv) <= 2:
        print("Give the names of two files with data. One at the begging of" +
              " simulation and the second at the end.")
        sys.exit()
    try:
        l1 = re.split(" ", ln.getline("InputParameters.csv", 9))
        spp_num = int(l1[2])
        l2 = re.split(" ", ln.getline("InputParameters.csv", 6))
        bitfit = int(l2[2])
        print("No. of pathogen species = " + str(spp_num) + "\n" +
              " Length of the bit string =" +
              str(int(re.split(" ", ln.getline("InputParameters.csv", 5))[2])))
        print("Length of bit string fit = " + str(bitfit))
    except Exception:
        print("Can't load the param file! You may be in a wrong directory.")
        sys.exit()
    try:
        L_init = asim.loadThePopulationBitstrings(sys.argv[1])
        print("First file loaded!")
    except Exception:
        print("Can't load file named " + str(sys.argv[1]) +
              ". Check if it exists.")
        sys.exit()
    try:
        L_endd = asim.loadThePopulationBitstrings(sys.argv[2])
        cloneList = asim.loadIndvPathoTags(sys.argv[2])
        print("Second file loaded!")
    except Exception:
        print("Can't load file named" + str(sys.argv[2]) +
              ". Check if it exists.")
        sys.exit()
    F_init = asim.hamDistInterSpecies(L_init)
    print("Similarities in the First file have been calculated!")
    F_endd = asim.hamDistInterSpecies(L_endd)
    print("Similarities in the Second file have been calculated!")
    clonez = asim.countClonesInSpeciesFromTags(cloneList[0])
    np.savetxt("cloneFreq.csv", clonez, fmt='%i')
    asim.plotCloneCount(clonez["numbOfIndv"], 25, len(cloneList[0]))
    print("Clonal variability of first pathogen population calculated!")
    # === More generic plot ===
    ax_label = 23
    T_label = 27
    TicksFS = 21
    transs = 0.8
    plt.figure(1, figsize=(24, 20))
    plt.subplot(221)
    xx = np.zeros(len(F_init))
    for ii, itm in enumerate(F_init):
        xx[ii] = np.mean(itm)
    plt.hist(xx, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("Inter-species similarity at the start of simulation",
              fontsize=T_label)
#    plt.xlabel("Inter-species similarity measure", fontsize=ax_label)
    plt.ylabel("Frequency of occurrence", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
#    plt.xlim(0., 1.)
    plt.subplot(222)
    xx = np.zeros(len(F_endd))
    for ii, itm in enumerate(F_endd):
        xx[ii] = np.mean(itm)
    plt.hist(xx, color=(0.3, 0.3, 0.3, transs), edgecolor="none")
    plt.title("Inter-species similarity at the end of simulation",
              fontsize=T_label)
#    plt.xlabel("Inter-species similarity measure", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
#    plt.xlim(0., 1.)
#    plt.savefig("SPP_sim_one.png")

    #  === Now the detailed plot! ===
    TicksFS_2 = TicksFS
    for ii in range(4, 7):
        plt.subplot(2, 3, ii)
        plt.hist(asim.hamDistInsideSpec(L_endd[ii], 15000),
                 color=(0.3, 0.3, 0.3, transs),  edgecolor="none")
        if (ii == 5):
            plt.title("Within-species similarity at the end of simulation"
                      + " (3 randomly selected species)", fontsize=T_label)
            plt.xlabel("Hamming distance", fontsize=ax_label)
        if (ii == 4):
            plt.ylabel("Frequency of occurrence", fontsize=ax_label)
        plt.xticks(fontsize=TicksFS_2)
        plt.yticks(fontsize=TicksFS_2)
        plt.grid(True)
    plt.savefig("SPP_sin_within_stop.png")


#    plt.show()

    print("DONE!")


if __name__ == "__main__":
    main()
