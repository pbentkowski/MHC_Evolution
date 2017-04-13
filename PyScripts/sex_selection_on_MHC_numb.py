#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Your doc string here, please...


Created on Thu Apr 13 15:42:41 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import numpy as np
import matplotlib as plt


def loadTheParents(moth="NumberOfMhcInMother.csv",
                   fath="NumberOfMhcInFather.csv"):
    """ """
    try:
        mother = np.genfromtxt(moth)[1::, 1::]
    except:
        print("Failed to load mothers MHC numbers. Check if file exists.")
        return None
    try:
        father = np.genfromtxt(fath)[1::, 1::]
    except:
        print("Failed to load fathers MHC numbers. Check if file exists.")
        return None
#    time = np.genfromtxt(moth)[:, 0]
    return mother, father


def trimData(mother, father, low=2, up=100):
    """ """
    new_mothers = []
    new_fathers = []
    for i, itm in enumerate(mother):
        if(np.mean(itm) >= low and np.mean(itm) <= up):
            new_mothers.append(itm)
            new_fathers.append(father[i])
    return np.array(new_mothers), np.array(new_fathers)


def reshapeMatherFather(mother, father):
    """ """
    if (mother.shape == father.shape):
        mother = np.reshape(mother, mother.shape[0] * mother.shape[1])
        father = np.reshape(father, father.shape[0] * father.shape[1])
        return mother, father
    else:
        print("Mother and father arrays need to have same shapes. Aborded.")
        return None


def pickMotherSizeGroups(mother, father):
    """ """
    ww = list(range(int(np.min(mother)), int(np.max(mother) + 1)))
    bigOnes = []
    for iii in ww:
        s_list = []
        for i, itm in enumerate(mother):
            if itm == iii:
                s_list.append(father[i])
#                s_list.append(itm - father[i])
        bigOnes.append(np.array(s_list))
    return ww, bigOnes


mother, father = ssmn.trimData(Mom, Fath, 2.5, 3.5)
rMom, rDad = ssmn.reshapeMatherFather(mother, father)
slope, interc, rval, pval, stderr = linregress(rMom, rDad)

x = rMom + 0.1 * np.random.randn(len(rMom))
y = rDad + 0.1 * np.random.randn(len(rDad))
xx = np.linspace(0, 9, 15)
ss = slope * xx + interc
plt.plot(x, y, '.')
plot(xx, ss, 'r-')
plt.show()