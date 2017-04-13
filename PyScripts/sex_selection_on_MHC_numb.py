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
import matplotlib.pyplot as plt
from scipy.stats import linregress


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
    """Takes reshaped 1D mother and father arrays"""
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


def plotAndDoStats(Mom, Dad, low_copy, up_copy):
    """ """
    if low_copy >= up_copy:
        print("Lower boudry on MHC copy number has to me smaller than higher",
              "boundry on MHC copy number.")
        return None
    mother, father = trimData(Mom, Dad, low_copy, up_copy)
    rMom, rDad = reshapeMatherFather(mother, father)
    slope, interc, rval, pval, stderr = linregress(rMom, rDad)
    print("\nslope =", slope, "; R^2 =", rval*rval, "; p_val =", pval)
    print("There are", len(rMom), "breeding pairs")
    x = rMom + 0.1 * np.random.randn(len(rMom))
    y = rDad + 0.1 * np.random.randn(len(rDad))
    xx = np.linspace(1, int(np.max(rMom)), 15)
    ss = slope * xx + interc
    FS = 16
    plt.figure(1, figsize=(8, 7))
    plt.plot(x, y, '.')
    plt.plot(xx, ss, 'r-')
    plt.grid(True)
    plt.xlabel("MHC copy number in 'mothers' (selecting indv.)", fontsize=FS)
    plt.ylabel("MHC copy number in 'fathers' (selected indv.)", fontsize=FS)
    plt.xticks(size=FS-1)
    plt.yticks(size=FS-1)
    plt.xlim(xmin=0.5)
    plt.ylim(ymin=0.5)
    plt.show()
    return pickMotherSizeGroups(rMom, rDad)
