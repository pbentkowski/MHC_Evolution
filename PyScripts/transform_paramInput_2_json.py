#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Iterates through the old simulation results and replaces the old
InputParameters.csv file with a newer Json version of the parameter input file,
that is easier to compare with templates in multi-simulation comparisons.

Created on Mon Jan 14 21:11:56 2019
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import json
import os
import re

paramz = ['alpha_factor_for_the_host_fitness_function',
          'heterozygote_advantage', 'host_gene_deletion_probability',
          'host_gene_duplication_probability',
          'host_maximal_number_of_genes_in_chromosome', 'host_population_size',
          'mutation_probability_in_host', 'mutation_probability_in_pathogen',
          'number_of_bits_per_antigen', 'number_of_bits_per_gene',
          'number_of_genes_per_host_one_chromosome',
          'number_of_host_generations',
          'number_of_pathogen_generation_per_one_host_generation',
          'number_of_pathogen_species', 'number_of_sex_mates',
          'number_of_threads', 'pathogen_population_size',
          'point_mutation_in_host_is_used', 'run_start_date_and_time',
          'separated_species_genomes']


def transformInputFile(path="."):
    """Reads data from a single InputParameters.csv file and loads them into
    the *.json version of it writing it in the same directory."""
    fileIn = os.path.join(path, "InputParameters.csv")
    lines = [line.rstrip('\n')[1::] for line in open(fileIn)]
    ww = {}
    for ii in lines:
        for kk in paramz:
            if(re.search(kk, ii.lower())):
                try:
                    ww[kk] = int(ii.split(" = ")[1])
                except Exception:
                    try:
                        ww[kk] = float(ii.split(" = ")[1])
                    except Exception:
                        ww[kk] = str(ii.split(" = ")[1])
    ww['number_of_sex_mates'] = 0
    ww['number_of_threads'] = 0
    fileOut = os.path.join(path, "InputParameters.json")
    with open(fileOut, 'w') as outfile:
        json.dump(ww, outfile, sort_keys=True, indent=4)


def doTheSwap(dirr=os.getcwd()):
    """Iterates through directories, finds InputParameters.csv and performs the
    swap function `transformInputFile()`. """
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.csv')):
                transformInputFile(dirName)


def main():
    """ """
    doTheSwap()


if __name__ == "__main__":
    main()
