#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Searches for 'InputParameters.json' files, pulls out parameters from them and
renders them in one line which can be feed as input to the model's program.

Created on Fri Jun 19 20:04:08 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import os
import sys
import json
import datetime as dt
# import linecache as ln


def readDate(string):
    """Takes a string and tries to convert it into a date. String has to have
    the ISO yyyy-mm-dd format"""
    try:
        dd = string.split("-")
        theDate = dt.date(int(dd[0]), int(dd[1]), int(dd[2]))
        return theDate
    except Exception:
        print("ERROR in readDate(): Bad string format! It has to be ISO's",
              "yyyy-mm-dd format!")
        return None


def loadTheDateFromParamFile(filePar):
    """Takes InputParameters.json file and tries to figure out what day the run
    was started."""
    try:
        with open(filePar) as f:
            prms = json.load(f)
        ll = prms['run_start_date_and_time'].split(".")[0].split("-")
#        ll = re.split(" ", ln.getline(filePar, 2))[2].split(".")[0].split("-")
    except Exception:
        print("ERROR in loadTheDate(): Cannot load the Params file. Check if",
              "the path to the params file is correct as well as its name.")
        return None
    try:
        theDay = dt.date(int(ll[0]), int(ll[1]), int(ll[2]))
        return theDay
    except Exception:
        print("ERROR in loadTheDate(): Cannot convert data into the date",
              "format. Check if the data file has the right flavour.")
        return None


def printTheParams3(theStartDate, dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type."""
    nogphoch = 'number_of_genes_per_host_one_chromosome'
    nopgpohg = 'number_of_pathogen_generation_per_one_host_generation'
    hmnogich = 'host_maximal_number_of_genes_in_chromosome'
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.json') and
               loadTheDateFromParamFile(filepath) >= theStartDate):
                ST = ""
                with open(filepath, 'r') as f:
                    prms = json.load(f)
                if prms['separated_species_genomes'] == "YES":
                    sepSpGen = "10"
                elif prms['separated_species_genomes'] == "NO":
                    sepSpGen = "11"
                else:
                    sepSpGen = prms['separated_species_genomes']
                ST = str(prms['number_of_threads']) + " "\
                    + str(prms['number_of_bits_per_gene']) + " "\
                    + str(prms['number_of_bits_per_antigen']) + " "\
                    + str(prms['host_population_size']) + " "\
                    + str(prms['pathogen_population_size']) + " "\
                    + str(prms['number_of_pathogen_species']) + " "\
                    + str(prms[nogphoch]) + " "\
                    + str(prms[nopgpohg]) + " "\
                    + str(prms['number_of_host_generations']) + " "\
                    + str(prms['mutation_probability_in_host']) + " "\
                    + str(prms['mutation_probability_in_pathogen']) + " "\
                    + sepSpGen + " "\
                    + str(prms['host_gene_deletion_probability']) + " "\
                    + str(prms['host_gene_duplication_probability']) + " "\
                    + str(prms[hmnogich]) + " "\
                    + str(prms['number_of_sex_mates']) + " "\
                    + str(prms['alpha_factor_for_the_host_fitness_function'])
                print(ST)


def main():
    """ """
    if len(sys.argv) <= 1:
        print("Give a starting date. It has to be in yyyy-mm-dd format.")
        sys.exit()
    try:
        startDate = readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert the input do a date format.")
        sys.exit()
    if startDate:
        printTheParams3(startDate)
    else:
        print("Wrong date format.")
        sys.exit()


if __name__ == "__main__":
    main()
