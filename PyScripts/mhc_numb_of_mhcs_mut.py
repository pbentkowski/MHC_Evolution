#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Tue Apr 28 13:32:49 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import os
import re
import pylab as p
import linecache as ln


def LoadGeneNumbers(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'HostsGeneDivers.csv'):
            genes = p.genfromtxt(filepath)
            arg.append((genes[:, 0], genes[:, 3], genes[:, 4], genes[:, 5],
                        genes[:, 1]))


def LoadTags(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'InputParameters.csv'):
                l = re.split(" ", ln.getline(filepath, 15))
                theValue = float(l[2])
                l = re.split(" ", ln.getline(filepath, 16))
                hh = l[2].split()[0]
                arg.append((theValue, hh))


def main():
    GeneData = []
    os.path.walk(os.getcwd(), LoadGeneNumbers, GeneData)
    Tags = []
    os.path.walk(os.getcwd(), LoadTags, Tags)

    AxLabelFontSize = 22
    AxisTickFontSize = 22
    AnnotateFontSize = 19

    annotScale = 10
    annotShift = 1

    Xmax = 2000
    Ymax = 120

    dec_places = '%1.3f'

    p.figure(1, figsize=(14, 7))
    i = 1
    for item in GeneData:
        XX = float(item[0][annotShift + i*annotScale])
        YY = float(item[1][annotShift + i*annotScale])
#        print XX, YY
        if Tags[i-1][1] == "NO":
            ax = p.plot(item[0], item[1], 'r-')
        else:
            ax = p.plot(item[0], item[1], 'b-')
        ax = p.annotate(dec_places % (Tags[i-1][0],), xy=(XX, YY),
                        xycoords='data', fontsize=AnnotateFontSize)
        i = i + 1
    p.ylabel('number of alles', fontsize=AxLabelFontSize)
    p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
    p.axis([0, Xmax, 0, Ymax])
    p.xticks(size=AxisTickFontSize)
    p.yticks(size=AxisTickFontSize)
    p.grid()

    p.figure(2, figsize=(14, 7))
    i = 0
    for item in GeneData:
        XX = float(item[0][annotShift + i*annotScale])
        YY = float(item[2][annotShift + i*annotScale])
        if Tags[i][1] == "NO":
            ax = p.plot(item[0], item[2], 'r-')
        else:
            ax = p.plot(item[0], item[2], 'b-')
        ax = p.annotate(dec_places % (Tags[i][0],), xy=(XX, YY),  xycoords='data',
                        fontsize=AnnotateFontSize)
        i = i + 1
    p.ylabel('Shannon\'s index', fontsize=AxLabelFontSize)
    p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
    p.axis([0, Xmax, 0, 10])
    p.xticks(size=AxisTickFontSize)
    p.yticks(size=AxisTickFontSize)
    p.grid()

    p.figure(3, figsize=(14, 7))
    i = 1
    for item in GeneData:
        XX = float(item[0][annotShift + i*annotScale])
        YY = float(item[3][annotShift + i*annotScale]/item[4][0])
#        print XX, YY
        if Tags[i-1][1] == "NO":
            ax = p.plot(item[0], item[3]/item[4], 'r-')
        else:
            ax = p.plot(item[0], item[3]/item[4], 'b-')
        ax = p.annotate(dec_places % (Tags[i-1][0],), xy=(XX, YY),
                        xycoords='data', fontsize=AnnotateFontSize)
        i = i + 1
    p.ylabel("fraction of homozygotes", fontsize=AxLabelFontSize)
    p.xlabel('time [host generations]', fontsize=AxLabelFontSize)
    p.axis([0, Xmax, 0, 1.0])
    p.xticks(size=AxisTickFontSize)
    p.yticks(size=AxisTickFontSize)
    p.grid()

    p.show()


if __name__ == "__main__":
    main()
