#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Reads files with genome size histograms ("HostGeneNumbTotal_ChrOne.csv" and
"HostMHCsNumbUniq_ChrOne.csv") and a general host statistics file
"HostsGeneDivers.csv" to render an animation how the size of the hosts'
chromosomes and the number of MHC unique alleles in them evolve.
You probably need a video codec like i.g. ffmpeg for matplotlib to be able
to create a MP4 animation clip.

Created on Tue Jun 30 12:37:54 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# This variable sets that every Nth row is used for plotting
# thanks to that we get a speedy animation.
everyOtherRow = 10

MAXX = 110
scale = 10**5

FontSize = 20
TickSize = 17
Y_max = 1.0
Hpad = 0.05
plt.rc('xtick', labelsize=TickSize)
plt.rc('ytick', labelsize=TickSize)


genMeans = np.genfromtxt("HostGeneNumbTotal_ChrOne.csv")
mhcMeans = np.genfromtxt("HostMHCsNumbUniq_ChrOne.csv")
GenerData = np.genfromtxt("HostsGeneDivers.csv")
print "Done loading data files!"

# -- trimmig rows
lastRow = genMeans[-1, :]
genMeans = genMeans[::everyOtherRow]
genMeans = np.vstack([genMeans, lastRow])
lastRow = mhcMeans[-1, :]
mhcMeans = mhcMeans[::everyOtherRow]
mhcMeans = np.vstack([mhcMeans, lastRow])
lastRow = GenerData[-1, :]
GenerData = GenerData[::everyOtherRow]
GenerData = np.vstack([GenerData, lastRow])
# -- trimmed

binz = np.arange(0.5, MAXX+0.5, 1.0)
tick_binz = np.arange(0, MAXX+1, 10.0)

fig = plt.figure(1, figsize=(15, 9))
ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
ax1.axis([0, GenerData[:, 0].max(), 0, 300])
ax1.set_ylabel('number of \nMHC allels', fontsize=FontSize)
ax1.set_xlabel('time (host generations)', fontsize=FontSize)
ax1.plot(GenerData[:, 0], GenerData[:, 3], 'k-')
ax1.grid(True)

ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=2)
plt.tight_layout(h_pad=Hpad)

line1, = ax1.plot([], [], 'ro', ms=10)
line3, = ax2.plot([], [], linewidth=2, color='k')


def init():
    line1.set_data([])
    return line1


def animate(i):
    line1.set_data(GenerData[i, 0], GenerData[i, 3])
    ax2.cla()
    ax2.axis([-1.0, 1.0, 0, 1.0])
    ax2.set_xlabel('number of genes per chromosome', fontsize=FontSize)
    ax2.set_ylabel('number of occurrences', fontsize=FontSize)
    plt.hist(genMeans[i, 1::], bins=binz, color=(0.3, 0.3, 0.3, 1.0),
             edgecolor="none")
    plt.hist(mhcMeans[i, 1::], bins=binz, color=(0.8, 0.0, 0.0, 1.0),
             edgecolor="none")
    plt.vlines(100, 0, 900, color="b", lw=2)
    plt.ylim(ymax=900)
    plt.xlim(xmax=MAXX)
    plt.xticks(tick_binz)
    plt.grid(True)
    plt.tight_layout(h_pad=Hpad)
    return line1,

ani = animation.FuncAnimation(fig, animate, frames=len(GenerData),
                              interval=200, blit=False)
ani.save('MHC_evol.mp4', fps=30)
plt.show()
