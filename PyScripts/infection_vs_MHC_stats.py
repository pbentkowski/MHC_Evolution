#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Your doc string here, please...


Created on Fri Feb  3 15:30:05 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def the_line(x, a, b):
    return a * x + b


def loadHostPopulation(FILE):
    '''Takes the file with all the hosts data loads it to a list. Each
    individual is loaded as a list of bit strings.And the population is a list
    of individuals.'''
    LL = []
    ll = []
    nextHost = False
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    if nextHost:
                        LL.append(ll)
                        ll = []
                else:
                    nextHost = True
                    ll.append(line.split()[0])
        LL.append(ll)
#        print j + 1
        return LL
    except IOError as e:
        print("I/O error({0}) in".format(e.errno),
              "loadTheHostPopulation(): {0}".format(e.strerror))


def loadPathoExposed(FILE):
    """ """
    LL = []
    try:
        with open(FILE) as infile:
            for line in infile:
                if re.search(r"#", line):
                    continue
                elif re.search(r"===", line):
                    ll = line.split()
                    if int(ll[7]) > 0:
                        LL.append(ll[11:-1])
                    else:
                        LL.append([])
                else:
                    continue
        return LL
    except IOError as e:
        print("I/O error({0}) in".format(e.errno),
              "loadPathoExposed(): {0}".format(e.strerror))


def uniqueMhcInHostOnly(hostPopList):
    """ """
    uniqHosts = []
    for indv in hostPopList:
        ll = []
        for mhc in indv:
            if mhc not in ll:
                ll.append(mhc)
        uniqHosts.append(ll)
    return uniqHosts


def calculateTheNumbers(hostPopList, pathoExposed):
    """ """
    uniq = uniqueMhcInHostOnly(hostPopList)
    LL = []
    for indv in uniq:
        LL.append(len(indv))
    uniqNumb = np.array(LL)
    ll = []
    for ii in pathoExposed:
        ll.append(len(ii))
    pathoNumb = np.array(ll)
    return uniqNumb, pathoNumb


def plotMHCvsPathoPresent(uniqNumb, pathoNumb, jitter=0.05):
    """ """
    jitterX = jitter * np.random.randn(uniqNumb.shape[0])
    jitterY = jitter * np.random.randn(pathoNumb.shape[0])
    slope, intercept, r_val, p_val, std_err = linregress(uniqNumb, pathoNumb)
    print("R^2 =", r_val**2)
    print("p =", p_val)
    FS = 16
    TicksFS = 14
    plt.figure(1, figsize=(12, 8))
    plt.scatter(uniqNumb + jitterX, pathoNumb + jitterY)
    plt.plot(uniqNumb, the_line(uniqNumb, slope, intercept), 'r-')
    plt.xlabel("number of unique MHC alleles in individual", fontsize=FS)
    plt.ylabel("number of infections presented to immune system", fontsize=FS)
    plt.xticks(fontsize=TicksFS)
    plt.yticks(fontsize=TicksFS)
#    plt.grid(True)
    plt.show()
