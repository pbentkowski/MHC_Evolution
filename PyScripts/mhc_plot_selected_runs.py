#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Plots generic statistics referring to time evolution of MHCs in hosts. These
are the number of MHC types, diversity of MHC types, host fitness etc.

Created on Mon May 18 17:03:31 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import os
import sys
import json
import pylab as p
# import linecache as ln


def meanMhcTypeNumber(DATA):
    """Calculates the mean number of MHC types (alleles) in the run after the
    500th time step."""
    ii = 0.
    mm = 0.
    for item in DATA:
        ii += 1.
        mm += p.mean(item[4][500::])
    return mm / ii


def loadTheData3(DIRR=os.getcwd()):
    """Iterates trough directories and look for HostsGeneDivers.csv file and
    corresponding InputParameters.json then copies the necessary data to
    a Python list to be analysed and plotted later on in the program. Version
    for Python 3.5 utilazing os.walk() function"""
    TheData = []
    for dirName, subdirList, fileList in os.walk(DIRR):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if filepath == os.path.join(dirName, 'HostsGeneDivers.csv'):
                genes = p.genfromtxt(filepath)
                paramsFile = os.path.join(dirName, 'InputParameters.json')
                with open(paramsFile) as f:
                    prms = json.load(f)
                # change here
                interestOne = float(prms['mutation_probability_in_pathogen'])
                # change here
                interestTwo = float(prms['number_of_pathogen_species'])
                path_spp = prms['number_of_pathogen_species']
                pop_size = float(prms['host_population_size'])
                if path_spp == "NOT_IN_THIS_MODEL":
                    path_spp = 1
                print("patho species:", path_spp, "| things:",  interestOne,
                      " ; ", interestTwo, "| dir:", dirName.split("/")[-1])
                TheData.append((interestOne, interestTwo, int(path_spp),
                                genes[:, 0], genes[:, 3], genes[:, 4],
                                genes[:, 5], genes[:, 2], genes[:, 6],
                                pop_size))
    return TheData


def main():
    """ """
#    TheData = []
#    os.path.walk(os.getcwd(), LoadTheData, TheData)
#    print(os.getcwd())

    TheData = loadTheData3()
    print(TheData)
    AxLabelFontSize = 22
    AxisTickFontSize = 22
    AnnotateFontSize = 19

    annotScale = 10
    annotShift = 200

    Xmax = 3000
    Ymax = 140
    textXlocal = 1500
    try:
        interestOne = float(sys.argv[1])
        interestTwo = float(sys.argv[2])
    except Exception:
        print("Can't recognise the parameters. Two are required.")
        sys.exit()
    saveFiggs = True  # True to save figures to disk, False to not save

    nnn = "One thing: " + str(interestTwo) + " Two thing: " + str(interestOne)

    dec_places = '%1.0f'

    p.figure(1, figsize=(14, 7))
    i = 1
    mm = 0.
    ii = 0.
    for item in TheData:
        print(item[1], interestTwo, item[0], interestOne)
        if (item[1] == interestTwo and item[0] == interestOne):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[4][annotShift + i*annotScale])
            p.plot(item[3], item[4], 'r-')
            mm += p.mean(item[4][500::])
            ii += 1.
            i = i + 1
            p.ylabel('number of  MHCs alles', fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, Ymax])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
#    p.hlines(mm / ii, 500, Xmax, colors='k', linestyles='solid', lw=2)
    ax = p.annotate(nnn, xy=(textXlocal, 180), xycoords='data',
                    fontsize=AnnotateFontSize)
    p.grid(True)
#    p.tight_layout()
    if saveFiggs:
        p.savefig("one_" + str(interestOne) + ".two_" +
                  str(interestTwo) + "_allel_num.png")

    p.figure(2, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == interestTwo and item[0] == interestOne):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[5][annotShift + i*annotScale])
            if item[2] == "NO":
                ax = p.plot(item[3], item[5], 'r-')
            else:
                ax = p.plot(item[3], item[5], 'b-')
            ax = p.annotate(dec_places % (item[0],), xy=(XX, YY),
                            xycoords='data', fontsize=AnnotateFontSize)
            i = i + 1
            p.ylabel('Shannon\'s index', fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 9])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    ax = p.annotate(nnn, xy=(textXlocal, 3.5), xycoords='data',
                    fontsize=AnnotateFontSize)
    p.grid()
#    p.tight_layout()
    if saveFiggs:
        p.savefig("one_" + str(interestOne) + ".two_" +
                  str(interestTwo) + "_Shann.png")

    p.figure(3, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == interestTwo and item[0] == interestOne):
            if item[2] == "NO":
                ax = p.plot(item[3], item[8]/item[6], 'r-')
            else:
                ax = p.plot(item[3], item[8]/item[6], 'b-')
            i = i + 1
            p.ylabel("CV fitness", fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 3.5])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
#    ax = p.annotate(nnn, xy=(textXlocal, 2.0), xycoords='data',
#                    fontsize=AnnotateFontSize)
    p.grid()
#    p.tight_layout()
    if saveFiggs:
        p.savefig("one_" + str(interestOne) + ".two_" +
                  str(interestTwo) + "_H_CV_fitt.png")

    p.figure(4, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == interestTwo and item[0] == interestOne):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[6][annotShift + i*annotScale])
            if item[2] == "NO":
                ax = p.plot(item[3], item[8]/float(item[2]), 'r-')
            else:
                ax = p.plot(item[3], item[8]/float(item[2]), 'b-')
            i = i + 1
            p.ylabel("hosts fitness", fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 4.0])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    p.grid()
#    p.tight_layout()
    if saveFiggs:
        p.savefig("one_" + str(interestOne) + ".two_" +
                  str(interestTwo) + "_H_fitt.png")

    p.figure(5, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == interestTwo and item[0] == interestOne):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[5][annotShift + i*annotScale])
            if item[2] == "NO":
                ax = p.plot(item[3], item[7] / item[9], 'r-')
            else:
                ax = p.plot(item[3], item[7] / item[9], 'b-')
#            ax = p.annotate(dec_places % (item[0],), xy=(XX, YY),
#                            xycoords='data', fontsize=AnnotateFontSize)
            i = i + 1
            p.ylabel('number of MHC loci / indv.', fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 100])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    ax = p.annotate(nnn, xy=(textXlocal, 3.5), xycoords='data',
                    fontsize=AnnotateFontSize)
    p.grid()
#    p.tight_layout()
    if saveFiggs:
        p.savefig("one_" + str(interestOne) + ".two_" +
                  str(interestTwo) + "AvgGenomeSize.png")

    p.show()


if __name__ == "__main__":
    main()
