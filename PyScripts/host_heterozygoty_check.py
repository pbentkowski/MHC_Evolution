# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 01:50:53 2016

@author: piotr
"""
import re
import sys


def loadHostPopulation(FILE):
    '''Takes the file with all the hosts data loads it to a list. Each
    individual is loaded as a list of gene tags. And the population is a list
    of individuals.'''
    B_list = []
    Gene_list = []
    contin = False
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif (re.search(r"===", line) and contin):
                    Gene_list.append(B_list)
                    B_list = []
                elif (re.search(r"===", line) and contin == False):
                    continue
                else:
                    contin = True
                    LL = line.split()
                    B_list.append(int(LL[3]))
        Gene_list.append(B_list)
        return Gene_list
    except IOError as e:
        print "I/O error({0}) in".format(e.errno),
        print "loadTheHostPopulation(): {0}".format(e.strerror)


def checkHeteroZyg(HH):
    """Checks how many of individuals contain only unique MHC alleles.
    Returns a fraction of population."""
    ii = 0
    for itm in HH:
        if len(itm)==len(set(itm)):
            ii += 1
    return float(ii) / float(len(HH))


def main():
    """The main function."""
    hh = loadHostPopulation(sys.argv[1])
    print "No reapets in", checkHeteroZyg(hh), "of genomes."


if __name__ == "__main__":
    main()
