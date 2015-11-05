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

open("dummy.5000.csv", 'w').close()
with open("PathoGenomesFile.3000.csv") as infile:
    for line in infile:
        ff = open("dummy.3000.csv", 'a')
        if re.search(r"#", line):
            ff.write(line)
        elif re.search(r"===", line):
            ff.write(line)
        else:
            ff.write(line.split()[0] + '\n')
        ff.close()
