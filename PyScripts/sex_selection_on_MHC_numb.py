#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Your doc string here...


Created on Thu Apr 13 15:42:41 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def loadTheParents(moth="NumberOfMhcInMother.csv",
                   fath="NumberOfMhcInFather.csv",
                   beforeMating="NumberOfMhcBeforeMating.csv"):
    """Simply load the data into two Numpy arrays."""
    try:
        mother = np.genfromtxt(moth)[1::, 1::]
    except Exception:
        print("Failed to load mothers MHC numbers. Check if file exists.")
        return None
    try:
        father = np.genfromtxt(fath)[1::, 1::]
    except Exception:
        print("Failed to load fathers MHC numbers. Check if file exists.")
        return None
    try:
        mates = np.genfromtxt(beforeMating)[1::, 1::]
    except Exception:
        print("Failed to load fathers MHC numbers. Check if file exists.")
        return None
#    time = np.genfromtxt(moth)[:, 0]
    return mother, father, mates


def trimData(mother, father, mates, low=0, up=100):
    """Takes data loaded by `loadTheParents()` and excludes time steps where
    average number of unique MHC types in mothers is lower or higher then
    user-defined limit."""
    new_mothers = []
    new_fathers = []
    new_mates = []
    for i, itm in enumerate(mother):
        if(np.mean(itm) >= low and np.mean(itm) <= up):
            new_mothers.append(itm)
            new_fathers.append(father[i])
            new_mates.append(mates[i])
    return np.array(new_mothers), np.array(new_fathers), np.array(new_mates)


def avrgMateMHCnumb(mate):
    """For each time step of the `mate` array it calculates the mean number
    of MHC types per time step and files an array of the exact shape as the the
    `mate` array. Used later for calculation."""
    mmFarh = np.zeros(mate.shape)
    for i, itm in enumerate(mate):
        mmFarh[i, :] = np.mean(itm)
    return mmFarh


def reshapeMatherFather(mother, father, mmFarh):
    """Simply reshapes the the data, but keep corresponding pairs."""
    if (mother.shape == father.shape and father.shape == mmFarh.shape):
        mother = np.reshape(mother, mother.shape[0] * mother.shape[1])
        father = np.reshape(father, father.shape[0] * father.shape[1])
        mmFarh = np.reshape(mmFarh, mmFarh.shape[0] * mmFarh.shape[1])
        return mother, father, mmFarh
    else:
        print("Mother and father arrays need to have same shapes. Aborded.")
        return None


def pickMotherSizeGroups(motherR, fatherR, mmFarhR):
    """Takes reshaped 1D mother and father arrays done by
    `reshapeMatherFather()`. For 'mothers' with N MHC types (`ww` list) it
    creates an array of all values of MHC types numbers of 'fathers'."""
    ww = list(range(int(np.min(motherR)), int(np.max(motherR) + 1)))
    bigOnes = []
    meanOnes = []
    for iii in ww:
        s_list = []
        m_list = []
        for i, itm in enumerate(motherR):
            if itm == iii:
                s_list.append(fatherR[i])
                m_list.append(mmFarhR[i])
#                s_list.append(itm - father[i])
        bigOnes.append(np.array(s_list))
        meanOnes.append(np.array(m_list))
    return ww, bigOnes, meanOnes


def meanFatherMHCnumb(ww, bigOnes, meanOnes):
    """ """
    meanFathr = np.zeros((len(ww), 3))
    for ii, moth in enumerate(ww):
        meanFathr[ii, 0] = moth
        meanFathr[ii, 1] = np.mean(bigOnes[ii])
        meanFathr[ii, 2] = np.mean(meanOnes[ii])
    return meanFathr


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
    plt.figure(1, figsize=(12, 11))
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


def plotDeviantFromMeanFather(mother, father, mate, lower, upper):
    """ """
    mother, father, mate = trimData(mother, father, mate, lower, upper)
    mmMate = avrgMateMHCnumb(mate)
    rMom, rDad, rMmMate = reshapeMatherFather(mother, father, mmMate)
    ww, bigOnes, meanOnes = pickMotherSizeGroups(rMom, rDad, rMmMate)
#    meanFathr = meanFatherMHCnumb(ww, bigOnes, meanOnes)
    ll = []
    for i, it in enumerate(ww):
        ll.append((it, np.mean(bigOnes[i] - meanOnes[i])))
    xx = np.array(ll)
    msSize = np.zeros(len(bigOnes))
    for i, itm in enumerate(bigOnes):
        msSize[i] = np.sqrt(float(len(itm)))
    plt.figure(1, figsize=(9, 6))
#    plt.plot(xx[:, 0], xx[:, 1], "o", ms=msSize)
    plt.scatter(xx[:, 0], xx[:, 1], s=msSize)
    plt.plot(xx[:, 0], np.zeros(len(xx)), 'k-', lw=2)
    plt.grid(axis='y')
    plt.xlabel("Number of MHC types in 'mothers'")
    plt.ylabel("Average deviation of 'fathers' MHC type number from" +
               " pre-mating population")
    plt.tight_layout()
    plt.savefig("SexSelectStrght.png")
    plt.show()
