#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 18:08:10 2015

@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
sys.setrecursionlimit(20000)


def probOfStreak(numBits, minFits, succProb=0.5, saved=None):
    """ """
    if saved is None:
        saved = {}
    ID = (numBits, minFits, succProb)
    if ID in saved:
        return saved[ID]
    else:
        if minFits > numBits or numBits <= 0:
            result = 0
        else:
            result = succProb**minFits
            for firsTrail in xrange(1, minFits+1):
                pr = probOfStreak(numBits-firsTrail, minFits, succProb, saved)
                result += (succProb**(firsTrail-1))*(1-succProb)*pr
        saved[ID] = result
    return result


def main():
    """ """
    try:
        numBits = int(raw_input("Enter length of the big bit string: "))
    except:
        print "The length of the bit string has to an integer!"
        sys.exit()
    try:
        minFits = int(raw_input("Enter length of the fitting string: "))
    except:
        print "The length of the fitting string has to an integer!"
        sys.exit()
    print "===================="
    print "Probability that a random bit string of length", minFits, "will",
    print "fit into a larger random bit string of length", numBits, "is:"
    print " p =", probOfStreak(numBits, minFits)

if __name__ == "__main__":
    main()
