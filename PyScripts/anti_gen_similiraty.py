#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Fri Oct 23 17:05:28 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import re
from bitstring import BitArray


def loadThePopulation(FILE):
    ''' '''
    LL = []
    spp_flag = False
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"===", line):
                    ll = []
                    spp_flag = True
                    print line,
                elif spp_flag:

    except IOError as e:
        print "I/O error({0}) in loadThePopulation(): {1}".format(e.errno,
                                                                  e.strerror)


def main():
    if len(sys.argv) <= 1:
        print "Give the name of the file with data."
        sys.exit()
    loadThePopulation(sys.argv[1])
    print "DONE!"


if __name__ == "__main__":
    main()
