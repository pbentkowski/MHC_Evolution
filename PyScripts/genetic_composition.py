#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Wed May 20 17:00:37 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import sys
import numpy as np
import matplotlib.pyplot as plt


def LoadParasites(File):
    genes = []
    counter = 0
    try:
        with open(File) as infile:
            for line in infile:
                if re.search(r"^(0|1)", line):
                    genes.append(int(line.split()[0], 2))
                    counter += 1
        return np.array(genes), counter
    except:
        print "ERROR in LoadParasite(): Can't load the genes."
        return None, None


def LoadHosts(File):
    genes = []
    counter = 0
    try:
        with open(File) as infile:
            for line in infile:
                if re.search(r"^(0|1)", line):
                    genes.append(int(line.split(" | ")[0], 2))
                    genes.append(int(line.split(" | ")[1], 2))
                    counter += 1
        return np.array(genes), counter
    except:
        print "ERROR in LoadHosts(): Can't load the genes."
        return None, None


def main():
    try:
        P, cp = LoadParasites("PathoGenomesFile.3000.csv")
        H, ch = LoadHosts("HostGenomesFile.3000.csv")
    except:
        print "Could not load the data files."
        sys.exit()

    plt.figure(1)
    plt.subplot(211)
    plt.hist(P, bins=cp)
    plt.grid(True)

    plt.subplot(212)
    plt.hist(H, bins=ch)
    plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
