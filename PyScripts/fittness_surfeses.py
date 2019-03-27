#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates 3D plots of the fitness function:
  f = P * exp( -(alpha * #MHC)^2 )
Used to evaluate the fitness of the hosts depending on the number of presented
pathogens species (P), the number of unique MHC variants (#MHC) and penalty
factor (alpha).

Created on Wed Mar 13 20:57:34 2019

@author: Piotr Bentkowski :: bentkowski.piotr@gmail.com

"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

figNumb = 0


def plotFinttSurf(alpha, mhc_max=30, patho_max=70):
    """The 3D plot."""
    global figNumb
    figNumb += 1
    fig = plt.figure(figNumb, figsize=(10, 10))
    ax = fig.gca(projection='3d')
    patho = np.arange(0, patho_max, 1.0)
    mhc = np.arange(1, mhc_max, 1.0)
    PP, MHC = np.meshgrid(patho, mhc, sparse=True)
    Z = PP * np.exp(-(alpha * MHC)**2)
    ax.plot_surface(PP, MHC, Z, cmap=cm.coolwarm, linewidth=0,
                    antialiased=False)
    ax.set_xlabel(r'Number of presented pathogen species ($P$)')
    ax.set_ylabel(r'Number of unique MHC variants ($N$)')
    ax.set_zlabel('fitness')
    ax.set_title(r"Surface of the fitness function " +
                 r"$f = P  e^{-(\alpha  N)^2}$ with $\alpha$ = "
                 + str(alpha))
    plt.tight_layout()
    plt.show()


def main():
    """ """
    plotFinttSurf(0.02)
    plotFinttSurf(0.08)


if __name__ == "__main__":
    main()
