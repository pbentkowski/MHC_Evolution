#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Your doc string here...


Created on Thu Apr 13 15:42:41 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import packed_plots_of_MHC_alleles as ppma


def loadTheParents(moth="NumberOfMhcInMother.csv",
                   fath="NumberOfMhcInFather.csv",
                   beforeMating="NumberOfMhcBeforeMating.csv"):
    """Simply load the data into two Numpy arrays."""
    try:
        mother = np.genfromtxt(moth)[1::, 1::]
    except Exception:
        print("Failed to load mothers MHC numbers. Check if file exists.")
        return None, None, None
    try:
        father = np.genfromtxt(fath)[1::, 1::]
    except Exception:
        print("Failed to load fathers MHC numbers. Check if file exists.")
        return None, None, None
    try:
        mates = np.genfromtxt(beforeMating)[1::, 1::]
    except Exception:
        print("Failed to load available mates MHC numbers.",
              "Check if file exists.")
        return None, None, None
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
    mmMate = np.zeros(mate.shape)
    for i, itm in enumerate(mate):
        mmMate[i, :] = np.mean(itm)
    return mmMate


def reshapeMatherFather(mother, father, mmMate):
    """Simply reshapes the the data, but keep corresponding pairs."""
    if (mother.shape == father.shape and father.shape == mmMate.shape):
        mother = np.reshape(mother, mother.shape[0] * mother.shape[1])
        father = np.reshape(father, father.shape[0] * father.shape[1])
        mmMate = np.reshape(mmMate, mmMate.shape[0] * mmMate.shape[1])
        return mother, father, mmMate
    else:
        print("Mother, father and mates arrays need to have same shapes.",
              "Aborded.")
        return None


def pickMotherSizeGroups(motherR, fatherR, mmFarhR):
    """Takes reshaped 1D mother, father and mate arrays done by
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
    """Plot regression plot between number of MHC types a 'mother' has and the
    number 'father' has. Not very useful though :-/ """
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
    """Plots the difference between the number of MHC types 'father' had and
    the mean number of MHC types that individuals had in host population before
    mating (the mean for the pool of available mates) for each size of 'mother'
    MHC repertoire."""
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


def justPlotDeviantFromMeanFather(ww, deltas, bSize, path):
    """Does the same as `plotDeviantFromMeanFather()` only it does not
    calculate the stats on it self."""
    bSize = np.sqrt(bSize)  # Create marker list
    plt.figure(1, figsize=(9, 6))
    plt.scatter(ww, deltas, s=bSize)
    plt.plot(ww, np.zeros(len(ww)), 'k-', lw=2)
    plt.grid(axis='y')
    plt.xlabel("Number of MHC types in 'mothers'")
    plt.ylabel("Average deviation of 'fathers' MHC type\nnumber from" +
               " pre-mating population")
    plt.tight_layout()
    fileName = os.path.join(path, "SexSelectStrght.png")
    plt.savefig(fileName)
    plt.cla()


def getTheData(theStartDate, templateList, dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type."""
    datOut = []
    for dirName, subdirList, fileList in os.walk(dirr):
        for file in fileList:
            filepath = os.path.join(dirName, file)
            if(filepath == os.path.join(dirName, 'InputParameters.json') and
               ppma.loadTheDateFromParamFile(filepath) >= theStartDate):
                paramzList = ppma.loadParamSettings(filepath)
                if ppma.compareParams(templateList, paramzList):
                    print("Processing dir:", dirName, end=" ")
                    moPth = os.path.join(dirName, 'NumberOfMhcInMother.csv')
                    faPth = os.path.join(dirName, 'NumberOfMhcInFather.csv')
                    mPth = os.path.join(dirName, 'NumberOfMhcBeforeMating.csv')
                    mother, father, mate = loadTheParents(moPth, faPth, mPth)
                    moth, fath, mate = trimData(mother, father, mate, 2, 100)
                    mmMt = avrgMateMHCnumb(mate)
                    rMom, rDad, rMmMt = reshapeMatherFather(moth, fath, mmMt)
                    ww, Fatrs, meanM = pickMotherSizeGroups(rMom, rDad, rMmMt)
                    bSize = np.zeros(len(Fatrs))
                    for i, itm in enumerate(Fatrs):
                        bSize[i] = len(itm)
                    deltas = []
                    for i, it in enumerate(ww):
                        deltas.append(np.mean(Fatrs[i] - meanM[i]))
                    justPlotDeviantFromMeanFather(ww, deltas, bSize, dirName)
                    try:
                        xx = np.transpose(np.vstack((ww, np.array(deltas),
                                                     bSize)))
                    except Exception:
                        print(" - failed to stack the data! Check if the",
                              "input file sizes (e.g. line numbers) are OK.")
                        continue
                    datOut.append(xx)
                    print(" - done.")
    return datOut
