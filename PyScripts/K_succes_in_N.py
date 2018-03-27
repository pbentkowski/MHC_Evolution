#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Calculates the theoretical probability of getting a run of K or more successes
(heads) in a row in N Bernoulli trials (coin flips)?

http://www.askamathematician.com/2010/07/q-whats-the-chance-of-getting-a-run-of-k-successes-in-n-bernoulli-trials-why-use-approximations-when-the-exact-answer-is-known/

Created on Fri Nov 13 18:08:10 2015

@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys


def probOfStreak(numBits, minFits, succProb=0.5, saved=None):
    """Recursive function calculating the probability that a define K-long
    bit string will fit into a N-long random bitstring."""
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
            for firsTrail in range(1, minFits+1):
                pr = probOfStreak(numBits-firsTrail, minFits, succProb, saved)
                result += (succProb**(firsTrail-1))*(1-succProb)*pr
        saved[ID] = result
    return result


def main():
    """ """
    try:
        numBits = int(input("Enter length of the big bit string: "))
    except Exception:
        print("The length of the bit string has to be an integer!")
        sys.exit()
    try:
        minFits = int(input("Enter length of the fitting string: "))
    except Exception:
        print("The length of the fitting string has to be an integer!")
        sys.exit()
    print("====================")
    print("Probability that a random bit string of length", minFits, "will " +
          "fit into a larger random bit string of length", numBits, "is:" +
          " p =", probOfStreak(numBits, minFits))
    print("\nDONE!")


if __name__ == "__main__":
    sys.setrecursionlimit(200000)
    main()
