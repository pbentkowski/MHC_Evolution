#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uses the files `NumberOfMhcInMother.csv`, `NumberOfMhcInFather.csv`,
`NumberOfMhcBeforeMating.csv` and `InputParameters.json` to analyse the
strength of selection preference on partners' MHC type number depending on the
sexual selection scenario used.


Created on Thu Apr 13 15:42:41 2017
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import os
import sys
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


def justPlotDeviantFromMeanFather(ww, deltas, bSize, path, suffix=""):
    """Does the same as `plotDeviantFromMeanFather()` only it does not
    calculate the stats on it self."""
    FS = 16
    bSize = np.sqrt(bSize)  # Create marker list
    plt.figure(1, figsize=(9, 6))
    plt.scatter(ww, deltas, s=bSize)
    plt.plot(ww, np.zeros(len(ww)), 'k-', lw=2)
    plt.grid(axis='y')
    plt.xlabel("Number of MHC types in 'mothers'", fontsize=FS)
    plt.ylabel("Average deviation of 'fathers' MHC type\nnumber from" +
               " pre-mating population", fontsize=FS)
    plt.xticks(size=FS-2)
    plt.yticks(size=FS-2)
    plt.tight_layout()
    if suffix:
        strr = "SexSelectStrght" + suffix + ".png"
        fileName = os.path.join(path, strr)
    else:
        fileName = os.path.join(path, "SexSelectStrght.png")
    plt.savefig(fileName)
    plt.cla()


def getTheData(theStartDate, templateList, dirr=os.getcwd()):
    """Walking the dir using Python 3.5. Variable theStartDate has to be
    a datetime.date() data type. Each item in the `datOut` structure is the
    result of computing one simulation."""
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


def avgDatOut(datOut):
    """Takes the data structure produced by `getTheData()` and averages the
    results collapsing them to one data structure ready to plot with
    `justPlotDeviantFromMeanFather()`"""
    minn = np.inf
    maxx = 0
    for itm in datOut:
        new_min = np.min(itm[:, 0])
        new_max = np.max(itm[:, 0])
        if new_min < minn:
            minn = new_min
        if new_max > maxx:
            maxx = new_max
    ww = np.arange(minn, maxx+1, 1)
    ll = np.zeros((len(ww), 4))
    ll[:, 0] = ww
    for itm in datOut:
        for ii, w in enumerate(ww):
            try:
                ll[ii, 1] += itm[ii][1]
                ll[ii, 2] += 1
                ll[ii, 3] += itm[ii][2]
            except Exception:
                continue
    out = np.zeros((len(ll), 3))
    out[:, 0] = ll[:, 0]
    out[:, 1] = ll[:, 1] / ll[:, 2]
    out[:, 2] = ll[:, 3]
    return out


def main():
    """Main function - the script's main body."""
    if len(sys.argv) <= 3:
        print("Two arguments are needed:")
        print("  1. Give a starting date. It has to be in yyyy-mm-dd format.")
        print("  2. Give the path to template file.")
        print("  3. Give the plot file suffix.")
        sys.exit()
    startDate = None
    try:
        startDate = ppma.readDate(sys.argv[1])
    except ValueError:
        print("Cannot convert argument #1 to a date format.")
        sys.exit()
    if startDate:
        try:
            template = ppma.loadParamSettings(sys.argv[2])
        except Exception:
            print("Cannot load the template file. Exiting.")
            sys.exit()
        try:
            theData = getTheData(startDate, template, ".")
#            print(theData)
        except Exception:
            print("Failed to process the data. Some serious issues arose.",
                  "Check if the cut-off host generation for calculating stats",
                  "is smaller than the total number of host generations.")
            sys.exit()
        if len(theData):
            np.save("sexSelectStrgt" + sys.argv[3], theData)
            out = avgDatOut(theData)
            justPlotDeviantFromMeanFather(out[:, 0], out[:, 1], out[:, 2],
                                          ".", sys.argv[3])
        else:
            print("No data files matching the criterions were found.",
                  "Specify your template file.")
            sys.exit()
    else:
        print("Wrong date format.")
        sys.exit()


if __name__ == "__main__":
    main()
