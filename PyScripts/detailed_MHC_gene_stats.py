#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Takes the files `HostGenomesFile.N.csv` and calculates how individual number of
variants are related with the number of all MHC genes by all individual genomes
and chromosomes division into chromosomes.

Created on Tue Oct 29 23:42:32 2019

@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import re
import numpy as np
import pandas as pd

stats_dt = np.dtype([('chr_1', np.int), ('chr_2', np.int), ('unq_1', np.int),
                     ('unq_2', np.int), ('tot', np.int), ('unq_tot', np.int)])


def loadHostPopulation(FILE="HostGenomesFile.0.csv"):
    '''Takes the file with all the hosts data loads it to a list. Each
    individual is loaded as as two lists (chromosome one and chromosome two)
    of gene tags. And the population is a list of individuals.'''
    B_list = [[], []]
    Gene_list = []
    contin = False
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif (re.search(r"===", line) and contin):
                    Gene_list.append(B_list)
                    B_list = [[], []]
                elif (re.search(r"===", line) and contin is False):
                    contin = True
                    continue
                else:
                    LL = line.split()
                    if LL[1] == 'ch_one':
                        B_list[0].append(int(LL[3]))
                    elif LL[1] == 'ch_two':
                        B_list[1].append(int(LL[3]))
        Gene_list.append(B_list)
        return Gene_list
    except IOError as e:
        print("I/O error({0}) in".format(e.errno) +
              " loadTheHostPopulation(): {0}".format(e.strerror))


def analiseGeneContent(Gene_list):
    """Takes what function `loadHostPopulation(...)` has produced and counts
    how many genes are there in individual (chromosome one, two, total,
    unique etc)."""
    stats = np.zeros(len(Gene_list), dtype=stats_dt)
    for i, indv in enumerate(Gene_list):
        stats['chr_1'][i] = len(indv[0])  # all genes on chromosome one
        stats['chr_2'][i] = len(indv[1])  # all genes on chromosome two
        stats['unq_1'][i] = len(set(indv[0]))  # unique MHCs on chromosome one
        stats['unq_2'][i] = len(set(indv[1]))  # unique MHCs on chromosome two
        stats['tot'][i] = len(indv[0] + indv[1])  # total number of MHC genes
        stats['unq_tot'][i] = len(set(indv[0] + indv[1]))  # total unique MHCs
    return pd.DataFrame(stats)
