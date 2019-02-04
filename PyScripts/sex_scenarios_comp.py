#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uses a post-processed file (you may need to edit it manually) and creates a
nice boxplot based on the data in that file. Check the `datype` data type to
see what kind of file you need. Script `packed_plots_of_MHC_alleles.py` may be
useful in creating this file. E.g. is `Integr_16_1e5`.


Created on Thu Apr  6 16:44:53 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


datype = np.dtype([('VAR', 'f8'), ('Mode', 'S12'), ('meanAllel', 'f8'),
                   ('stdAllel', 'f8'), ('slope', 'f8'), ('indvMean', 'f8'),
                   ('indvSTD', 'f8'), ('meanFitt', 'f8'), ('cvFitMean', 'f8'),
                   ('meanPato', 'f8'), ('stdPato', 'f8')])


def plotBoxesGeneMeans(dataArr, suffix="", ymaxes=()):
    """Plots boxplots of how many MHC types are there in different simulation,
    and fitnesses from data given in files usualy called `DataSlice.csv`"""
    if len(ymaxes) < 3:
        ymaxes = ()
        print("Setting y-maxes to default")
    if dataArr.dtype == datype:
        lebls = np.unique(dataArr["Mode"])
    else:
        print("ERROR in plotHistograms(): wrong numpy data type. It should",
              "be:", datype)
        return None
    if suffix == "":
        figureNameOne = "sexComprMHC.png"
    else:
        figureNameOne = "sexComprMHC_" + suffix + ".png"
    lbls = []
    for itm in lebls:
        lbls.append(itm.decode())
    ll1 = []
    ll2 = []
    ll3 = []
    fs = 16
    tkfs = 16
    # MHC numbers in population and in individuals
    plt.figure(1, figsize=(18, 6))
    plt.subplot(131)
    for itm in lebls:
        ll1.append(dataArr[dataArr["Mode"] == itm]['meanAllel'])
    for itm in lebls:
        ll2.append(dataArr[dataArr["Mode"] == itm]['indvMean'])
    for itm in lebls:
        ll3.append(dataArr[dataArr["Mode"] == itm]['meanPato'])
    boxprops = dict(linestyle='-', linewidth=2.5, color='k')
    medianprops = dict(linestyle='-', linewidth=2.5)
    whiskerprops = dict(linewidth=2.5)
    capprops = dict(linewidth=2.5)
    flierprops = dict(markersize=10)
    plt.boxplot(ll1, labels=lbls, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops)
    plt.ylabel("number of types MHCs in population", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.ylim(bottom=0)
    if ymaxes:
        plt.ylim(top=ymaxes[0])
    plt.grid(axis='y')
#    plt.figure(2, figsize=(10, 8))
    plt.subplot(132)
    plt.boxplot(ll2, labels=lbls, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops)
    plt.ylabel("number of MHC in individual", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.ylim(bottom=0)
    if ymaxes:
        plt.ylim(top=ymaxes[1])
    plt.grid(axis='y')

    # mean normalised presented pathogens
    plt.subplot(133)
    plt.boxplot(ll3, labels=lbls, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops)
    plt.ylabel("mean normalised number of presented pathogens", fontsize=fs)
    plt.xticks(fontsize=tkfs)
    plt.yticks(fontsize=tkfs)
    plt.ylim(bottom=0)
    if ymaxes:
        plt.ylim(top=ymaxes[2])
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(figureNameOne)
    plt.show()


def main():
    """"Main function - the script's main body."""
    if len(sys.argv) <= 1:
        print("One argument is needed:")
        print("  1.  Give the name of the data file.")
        print(" [2]. Suffix added to the figure file.")
        sys.exit()
    try:
        dat = np.genfromtxt(sys.argv[1], dtype=datype)
#        print(dat.dtype.descr)
    except Exception:
        print("Cannot load the data from data file. Failed.")
        sys.exit()
#    if True:
    try:
        if len(sys.argv) > 2:
            maxes = (120, 70, 0.4)  # set these for plot axes max-mins
            # maxes = (populMHCnumb, IndvMHCnumb, PresentPatho)
            plotBoxesGeneMeans(dat, sys.argv[2], maxes)
            # sys.argv[2] is suffix added to the name of the PNG plot with data
        else:
            plotBoxesGeneMeans(dat)
#    else:
    except Exception:
        print("Failed to plot the data.")
        sys.exit()


if __name__ == "__main__":
    main()
