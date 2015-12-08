#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Thu Oct 29 17:47:26 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import sys

try:
    open(sys.argv[2], 'w').close()
except:
    print "Cannot create the new file"
    sys.exit()
try:
    with open(sys.argv[1]) as infile:
        for line in infile:
            ff = open(sys.argv[2], 'a')
            if re.search(r"#", line):
                ff.write(line)
            elif re.search(r"===", line):
                ff.write(line)
            else:
                ff.write(line.split()[0] + '\n')
            ff.close()
except:
    print "Cannot process the input file. Check if it exists."
    sys.exit()
