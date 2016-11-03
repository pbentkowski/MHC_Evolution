#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Your doc string here, please...


Created on Thu Nov  3 16:26:40 2016
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import numpy as np
import bitstring as bts


def generateBistring(length, pp=0.5):
    """ """
    bb = ""
    for ii in range(length):
        if np.random.rand() <= pp:
            bb += "1"
        else:
            bb += "0"
    bbit = bts.BitString(bin=bb)
    return bbit


def dumpLongBitToNumbers(longBit, frameSize):
    """ """
    sizze = len(longBit)-frameSize + 1
    arr = np.zeros(sizze, dtype=int)
    for ii in range(sizze):
        arr[ii] = longBit[ii:ii+frameSize].int
    return arr


def createMHCsArr(howMuch, length):
    """ """
    arr = np.zeros(howMuch, dtype=int)
    for ii in range(howMuch):
        arr[ii] = generateBistring(length).int
    return arr


def createAtingArr(howMuch, length, frameSize):
    """ """
    antiArr = []
    for ii in range(howMuch):
        antiArr.append(dumpLongBitToNumbers(generateBistring(length),
                                            frameSize))
    return antiArr


def main():
    """ """
    if len(sys.argv) <= 3:
        print("The script needs 3 arguments:")
        antiSize = int(input("  1. Size of the antigen (big bitstrig): "))
        mhcSize = int(input("  2. Size of the MHC (small bitstrig): "))
        theNumber = int(input("  3. Number of MHC to test: "))
    else:
        try:
            antiSize = int(sys.argv[1])
            mhcSize = int(sys.argv[2])
            theNumber = int(sys.argv[3])
            if(antiSize < mhcSize):
                print("Size of antigen cannot be smaller then size of MHC.",
                      "Quit.")
                sys.exit()
        except:
            print("Cannot convert arguments to numbers. Quit")
            sys.exit()

    Antig = createAtingArr(theNumber, antiSize, mhcSize)
    MHCs = createMHCsArr(theNumber, mhcSize)
    counter = 0.
    for itm in zip(MHCs, Antig):
        for ant in itm[1]:
            if itm[0] == ant:
                counter += 1.
                break
    print("MHC fitted", counter, "times")
    print("Fraction of MHC tested that fitted:",
          counter / theNumber)
    print("\nDone!")

if __name__ == "__main__":
    main()
