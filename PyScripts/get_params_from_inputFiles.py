#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Searches for 'InputParameters.csv' files, pulls out parameters from them and
renders them in one line which can be feed as input to the model's program.

Created on Fri Jun 19 20:04:08 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import os
#import re
#import pylab as p
#import linecache as ln


def LoadTheData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'InputParameters.csv'):
            strr = ""
            with open(filepath, 'r') as f:
                for ii, line in enumerate(f):
                    l = line.split()
                    if ii > 2 and len(l) > 1:
                        strr += l[2] + " "
            print strr


def main():
    TheData = []
    os.path.walk(os.getcwd(), LoadTheData, TheData)

if __name__ == "__main__":
    main()
