#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:27:47 2017

@author: piotr
"""
import re
import numpy as np
import matplotlib.pyplot as plt


def getAgeOfGenes(FILE, firstIndex=1):
    """ """
    times = []
    with open(FILE) as infile:
        i = 0
        for line in infile:
            if re.search(r"#", line):
                continue
            elif(i < firstIndex):
                continue
            else:
                LL = line.split()
                for j, itm in enumerate(LL):
                    if()

            i += 1
