#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Mon May 18 17:03:31 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import os
import re
import pylab as p
import linecache as ln


def LoadTheData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'HostsGeneDivers.csv'):
            genes = p.genfromtxt(filepath)
            paramsFile = os.path.join(dirname, 'InputParameters.csv')
            l = re.split(" ", ln.getline(paramsFile, 15))   # change here
            InterestigThing = float(l[2])
            l = re.split(" ", ln.getline(paramsFile, 11))
            geneNumb = float(l[2])
            l = re.split(" ", ln.getline(paramsFile, 16))
            hh = l[2].split()[0]
            print "antigens:", geneNumb, "| thing:", InterestigThing, " | ", hh
            arg.append((InterestigThing, geneNumb, hh, genes[:, 0],
                        genes[:, 3], genes[:, 4], genes[:, 5],
                        genes[:, 1], genes[:, 6]))


def main():
    TheData = []
    os.path.walk(os.getcwd(), LoadTheData, TheData)

    AxLabelFontSize = 22
    AxisTickFontSize = 22
    AnnotateFontSize = 19

    annotScale = 10
    annotShift = 200

    Xmax = 2000
    Ymax = 60

    pathoGenSize = 5  # change to select a different set of data
    pathoNumSpec = 0.001  # change to select a different set of data
    saveFiggs = False  # True to save figures to disk, False to not save

    nnn = "antigens: " + str(pathoGenSize) + "   muts: " + str(pathoNumSpec)

    dec_places = '%1.0f'

    p.figure(1, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == pathoGenSize and item[0] == pathoNumSpec):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[4][annotShift + i*annotScale])
    #        print XX, YY
            if item[2] == "NO":
                ax = p.plot(item[3], item[4], 'r-')
            else:
                ax = p.plot(item[3], item[4], 'b-')
#            ax = p.annotate(dec_places % (item[0],), xy=(XX, YY),
#                            xycoords='data', fontsize=AnnotateFontSize)
            i = i + 1
            p.ylabel('number of alles', fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, Ymax])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    ax = p.annotate(nnn, xy=(1400, 50), xycoords='data',
                    fontsize=AnnotateFontSize)
    p.grid()
    if saveFiggs:
        p.savefig("g_" + str(pathoGenSize) + ".s_" +
                  str(pathoNumSpec) + "_geneNum.png")

    p.figure(2, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == pathoGenSize and item[0] == pathoNumSpec):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[5][annotShift + i*annotScale])
    #        print XX, YY
            if item[2] == "NO":
                ax = p.plot(item[3], item[5], 'r-')
            else:
                ax = p.plot(item[3], item[5], 'b-')
            ax = p.annotate(dec_places % (item[0],), xy=(XX, YY),
                            xycoords='data', fontsize=AnnotateFontSize)
            i = i + 1
            p.ylabel('Shannon\'s index', fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 4])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    p.grid()
    if saveFiggs:
        p.savefig("g_" + str(pathoGenSize) + ".s_" +
                  str(pathoNumSpec) + "_Shann.png")

    p.figure(3, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == pathoGenSize and item[0] == pathoNumSpec):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[6][annotShift + i*annotScale]/item[7][0])
    #        print XX, YY
            if item[2] == "NO":
                ax = p.plot(item[3], item[6]/item[7], 'r-')
            else:
                ax = p.plot(item[3], item[6]/item[7], 'b-')
#            ax = p.annotate(dec_places % (item[0],), xy=(XX, YY),
#                            xycoords='data', fontsize=AnnotateFontSize)
            i = i + 1
            p.ylabel("fraction of homozygotes", fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 1.0])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    ax = p.annotate(nnn, xy=(200, 0.8), xycoords='data',
                    fontsize=AnnotateFontSize)
    p.grid()
    if saveFiggs:
        p.savefig("g_" + str(pathoGenSize) + ".s_" +
                  str(pathoNumSpec) + "_HZ.png")

    p.figure(4, figsize=(14, 7))
    i = 1
    for item in TheData:
        if (item[1] == pathoGenSize and item[0] == pathoNumSpec):
            XX = float(item[3][annotShift + i*annotScale])
            YY = float(item[8][annotShift + i*annotScale])
    #        print XX, YY
            if item[2] == "NO":
                ax = p.plot(item[3], item[8], 'r-')
            else:
                ax = p.plot(item[3], item[8], 'b-')
#            ax = p.annotate(dec_places % (item[0],), xy=(XX, YY),
#                            xycoords='data', fontsize=AnnotateFontSize)
            i = i + 1
            p.ylabel("hosts fitness", fontsize=AxLabelFontSize)
            p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
            p.axis([0, Xmax, 0, 50.0])
            p.xticks(size=AxisTickFontSize)
            p.yticks(size=AxisTickFontSize)
    ax = p.annotate(nnn, xy=(300, 45), xycoords='data',
                    fontsize=AnnotateFontSize)
    p.grid()
    if saveFiggs:
        p.savefig("g_" + str(pathoGenSize) + ".s_" +
                  str(pathoNumSpec) + "_CV.png")

    p.show()


if __name__ == "__main__":
    main()
