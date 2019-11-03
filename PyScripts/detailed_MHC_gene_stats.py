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
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import linregress
import packed_plots_of_MHC_alleles as ppma

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


def getTheGenes(theStartDate, templateList, dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type."""
    vv = ppma.lookForVARinList(templateList)
    datOut = []
#    dataOrdering = ['VAR', 'VARX', 'meanAllel', 'stdAllel', 'slope']
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.json') and
               ppma.loadTheDateFromParamFile(filepath) >= theStartDate):
                paramzList = ppma.loadParamSettings(filepath)
#                with open(filepath) as f:
#                    prms = json.load(f)
                if ppma.compareParams(templateList, paramzList):
                    print("Data from:", dirName, end=" ")
                    popFiles = os.path.join(dirName, "HostGenomesFile.*.csv")
                    for fil in glob.glob(popFiles):
                        if not re.search('HostGenomesFile.0.csv', fil):
                            hostPopFile = fil
                    Gene_list = loadHostPopulation(hostPopFile)
                    geneStats = analiseGeneContent(Gene_list)
                    var = float(paramzList[vv['VAR']])
                    varx = float(paramzList[vv['VARX']])
                    geneStats['spp'] = varx
                    geneStats['patho_mut'] = var
                    datOut.append((geneStats))
                    print("- done!")
    return datOut


def plotFraction(result):
    """ """
    result['uniqFrac'] = result['unq_tot'] / result['tot']
    TicksFS = 14
    xtks = np.unique(result['spp']).astype('int')
    ytks = np.arange(0, 11, 2) / 10
    fig1, ax1 = plt.subplots(1, 1, figsize=(12, 8))
    ax1 = sns.boxenplot(x="spp", y="uniqFrac", hue="patho_mut", data=result)
    ax1.legend_.remove()
    ax1.legend(title="Patho. mut.", title_fontsize=TicksFS-1,
               fontsize=TicksFS-2, loc=4, edgecolor='white')
    ax1.set_yticklabels(ytks, fontsize=TicksFS)
    ax1.set_ylabel("fraction of unique MHC variants in individual",
                   fontsize=TicksFS-1)
    ax1.set_xticklabels(xtks, fontsize=TicksFS)
    ax1.set_xlabel("number of pathogen species", fontsize=TicksFS-1)
    ax1.grid(True, axis="y")
    ax1.set_ylim((0, 1.05))
    ax1.set_title("Fraction of unique MHC variants in individual's genome",
                  fontsize=TicksFS+2)
    plt.tight_layout()
    fig1.savefig("fraction_of_unique_MHC.png")


def plotHetero(result):
    """ """
    result['hetero'] = result['unq_tot'] / (result['unq_1'] + result['unq_2'])
    TicksFS = 14
    xtks = np.unique(result['spp']).astype('int')
#    ytks = np.arange(0, 11, 2) / 10
    fig1, ax1 = plt.subplots(1, 1, figsize=(12, 8))
    ax1 = sns.boxenplot(x="spp", y="hetero", hue="patho_mut", data=result)
    ax1.legend_.remove()
    ax1.legend(title="Patho. mut.", title_fontsize=TicksFS-1,
               fontsize=TicksFS-2, loc=4, edgecolor='white')
#    ax1.set_yticklabels(ytks, fontsize=TicksFS)
    ax1.set_ylabel("heterozygosity approximation", fontsize=TicksFS-1)
    ax1.set_xticklabels(xtks, fontsize=TicksFS)
    ax1.set_xlabel("number of pathogen species", fontsize=TicksFS-1)
    ax1.grid(True, axis="y")
    ax1.set_ylim((0.45, 1.05))
    ax1.set_title("INV / (uniq_chr_one + uniq_chr_two) : range[0.5 , 1.0]",
                  fontsize=TicksFS+2)
    plt.tight_layout()
    fig1.savefig("hetero_apprx_MHC.png")


def plotTotNumb(result):
    """ """
    TicksFS = 14
    xtks = np.unique(result['spp']).astype('int')
#    ytks = np.arange(0, 11, 2) / 10
    fig3, ax3 = plt.subplots(1, 1, figsize=(12, 8))
    ax3 = sns.barplot(x="spp", y="tot", hue='patho_mut', data=result,
                      edgecolor='white', saturation=0.25)
    ax3 = sns.barplot(x="spp", y="unq_tot", hue='patho_mut', data=result,
                      edgecolor='white', saturation=0.8)
    ax3.legend_.remove()
    ax3.legend(title="Patho. mut.", title_fontsize=TicksFS-1,
               fontsize=TicksFS-2, loc=1, edgecolor='white')
#    ax3.set_yticklabels(ytks, fontsize=TicksFS)
    ax3.set_ylabel("number of MHCs in individual", fontsize=TicksFS-1)
    ax3.set_xticklabels(xtks, fontsize=TicksFS)
    ax3.set_xlabel("number of pathogen species", fontsize=TicksFS-1)
#    ax3.grid(True, axis="y")
    ax3.set_ylim((0, 10))
    ax3.set_title("MHCs in individual's genome (unique and all loci)",
                  fontsize=TicksFS+2)
    plt.tight_layout()
    fig3.savefig("number_MHC.png")


def chromVsUniq(result, chrom='chr_1', unq='unq_1'):
    """ """
    spp = np.unique(result['spp'])
    mut = np.unique(result['patho_mut'])
    nn = len(spp) * len(mut)
    out = np.zeros(nn, dtype=[('spp', np.int), ('patho_mut', np.float),
                              ('corr', np.float), ('p_val', np.float),
                              ('slope', np.float)])
    i = 0
    for mt in mut:
        for sp in spp:
            trim = result[(result['spp'] == sp) & (result['patho_mut'] == mt)]
            slp, inter, r_val, p_val, std_err = linregress(trim[chrom],
                                                           trim[unq])
#            print(slp, inter, r_val, p_val, std_err)
            out['spp'][i] = int(sp)
            out['patho_mut'][i] = mt
            out['corr'][i] = r_val
            out['p_val'][i] = p_val
            out['slope'][i] = slp
            i += 1
    return pd.DataFrame(out)


def plotChrVsUnq(result, what='corr', chrom='chr_1', unq='unq_1'):
    """ """
    chrVsUnq = chromVsUniq(result, chrom, unq)
    TicksFS = 14
    xtks = np.unique(result['spp']).astype('int')
#    ytks = np.arange(0, 11, 2) / 10
    fig4, ax4 = plt.subplots(1, 1, figsize=(12, 8))
    ax4 = sns.barplot(x="spp", y=what, hue='patho_mut', data=chrVsUnq)
    ax4.legend_.remove()
    ax4.legend(title="Patho. mut.", title_fontsize=TicksFS-1,
               fontsize=TicksFS-2, loc=1, edgecolor='white')
#    ax4.set_yticklabels(ytks, fontsize=TicksFS)
    ax4.set_ylabel("number of MHCs in individual", fontsize=TicksFS-1)
    ax4.set_xticklabels(xtks, fontsize=TicksFS)
    ax4.set_xlabel("number of pathogen species", fontsize=TicksFS-1)
#    ax4.grid(True, axis="y")
    if (what == 'corr'):
        ax4.set_ylabel("R correlation coef.", fontsize=TicksFS-1)
        ax4.set_ylim((0, 1.1))
        ax4.set_title("R correlation between unique variants and all MHC",
                      fontsize=TicksFS+2)
    if (what == 'slope'):
        ax4.set_ylabel("slope of correlation", fontsize=TicksFS-1)
        ax4.set_ylim((0, 1.5))
        ax4.set_title("Slope of correlation between unique variants"
                      + " and all MHC", fontsize=TicksFS+2)
    plt.tight_layout()
    figName = "Corr_" + str(what) + "_" + str(chrom) + "_" + str(unq) + ".png"
    fig4.savefig(figName)


def chromoProp(row):
    """ """
    minn = np.minimum(row['unq_1'], row['unq_2'])
    maxx = np.maximum(row['unq_1'], row['unq_2'])
    return minn / maxx


def plotChromoProp(result):
    """ """
    result['chro_prop'] = result.apply(lambda row: chromoProp(row), axis=1)
    TicksFS = 14
    xtks = np.unique(result['spp']).astype('int')
    ytks = np.arange(0, 11, 2) / 10
    fig1, ax1 = plt.subplots(1, 1, figsize=(12, 8))
    ax1 = sns.boxenplot(x="spp", y="chro_prop", hue="patho_mut", data=result)
    ax1.legend_.remove()
    ax1.legend(title="Patho. mut.", title_fontsize=TicksFS-1,
               fontsize=TicksFS-2, loc=4, edgecolor='white')
    ax1.set_yticklabels(ytks, fontsize=TicksFS)
    ax1.set_ylabel("proportion of the chromosome sizes",
                   fontsize=TicksFS-1)
    ax1.set_xticklabels(xtks, fontsize=TicksFS)
    ax1.set_xlabel("number of pathogen species", fontsize=TicksFS-1)
    ax1.grid(True, axis="y")
    ax1.set_ylim((0, 1.05))
    ax1.set_title("Proportion of sizes of chromosomes",
                  fontsize=TicksFS+2)
    plt.tight_layout()
    fig1.savefig("chrom_size_prop.png")


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 2:
        print("Two arguments are needed:")
        print("  1. Give a starting date. It has to be in yyyy-mm-dd format.")
        print("  2. Give the path to template file.")
        sys.exit()
    startDate = None
    try:
        startDate = ppma.readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert argument #1 to a date format.")
        sys.exit()
    if startDate:
        try:
            template = ppma.loadParamSettings(sys.argv[2])
            if template is None:
                print("Failed to load the template file. Exiting.",
                      "Check if the path is correct - you may wish to provide",
                      "an absolute path.")
                sys.exit()
        except Exception:
            print("Cannot load the template file. Exiting.")
            sys.exit()
        try:
            datOut = getTheGenes(startDate, template, os.getcwd())
            result = pd.concat(datOut, ignore_index=True)
            plotFraction(result)
            plotHetero(result)
            plotTotNumb(result)
            plotChrVsUnq(result, 'corr', 'chr_2', 'unq_2')
            plotChrVsUnq(result, 'slope', 'chr_2', 'unq_2')
            plotChrVsUnq(result, 'corr', 'chr_1', 'unq_1')
            plotChrVsUnq(result, 'slope', 'chr_1', 'unq_1')
            plotChromoProp(result)
        except Exception:
            print("Failed to process the data. Some serious issues arose.",
                  "Check if the cut-off host generation for calculating stats",
                  "is smaller than the total number of host generations.")
            sys.exit()


if __name__ == "__main__":
    main()
