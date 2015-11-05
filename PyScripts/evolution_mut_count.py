#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Thu Nov  5 14:30:59 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import sys
import linecache as ln
import numpy as np
import matplotlib.pyplot as plt
#import bitstring as bts


def loadHostPopulation(FILE):
    '''Takes the file with all the hosts data loads it to a list. Each
    individual is loaded as a list of bit strings.And the population is a list
    of individuals.'''
    B_list = []
    Mut_list = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    continue
                else:
                    LL = line.split()
                    bb = int(LL[3])
                    if bb in B_list:
                        pass
                    else:
                        B_list.append(bb)
                        Mut_list.append((len(LL) - 4) // 2)
        return np.array(Mut_list)
    except IOError as e:
        print "I/O error({0}) in".format(e.errno),
        print "loadTheHostPopulation(): {0}".format(e.strerror)


def main():
    if len(sys.argv) <= 1:
        print "Give the name of the file with data at the end of simulation."
        sys.exit()
    try:
        l = re.split(" ", ln.getline("InputParameters.csv", 9))
        print "No. of pathogen species =", int(l[2])
    except:
        print "Can't find parameter file! You may be in a wrong directory."
        sys.exit()
    try:
        L_endd = loadHostPopulation(sys.argv[1])
        print "File loaded!"
    except:
        print "Can't load file named", sys.argv[1], ". Check if it exists."
        sys.exit()
    # === More generic plot ===
    ax_label = 20
    T_label = 24
    TicksFS = 18
    transs = 0.8
    plt.figure(1, figsize=(12, 7))
    plt.hist(L_endd, color=(0.3, 0.3, 0.7, transs), edgecolor="none")
    plt.title("End of simulation", fontsize=T_label)
    plt.xlabel("Number of accumulated mutation in MHCs", fontsize=ax_label)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
    plt.grid(True)
    plt.xlim(0., 20.)
    plt.savefig("HOST_muts_end.png")

    plt.show()

    print "DONE!"


if __name__ == "__main__":
    main()
