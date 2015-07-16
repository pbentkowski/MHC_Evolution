#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Doc string here...

Created on Thu Mar 19 13:30:02 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznan, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

li = 0.65
lii = 0.75
FS = 16
TFS = 13


def secondMax(data):
    W = np.where(data[-1, 1:] == data[-1, 1:].max())[0][0] + 1
    data[:, W] = 0
    for i in np.arange(len(data))[::-1]:
        mx = data[i, 1:].max()
        if mx > 0:
            W1 = np.where(data[i, 1:] == mx)[0][0] + 1
            break
    return W1


def main():
    try:
        data = np.genfromtxt("PathoPopSizes.csv")
    except:
        print "No file with data or argumnets could be found!"
        sys.exit()
    # Find max:
    W = np.where(data[-1, 1:] == data[-1, 1:].max())[0][0] + 1
    MM = lambda x, M: M * np.ones(x.shape[0])
    mm = np.mean(data[0, 1::])

    fig1 = plt.figure(1, figsize=(14, 6))
    fig1.canvas.set_window_title('Parasite species populations sizes')
    for i in np.arange(1, data.shape[1]):
        plt.plot(data[:, 0], data[:, i], '-', lw=1, color=(li, li, li, 0.95))
        print "Spec. No", i, "is", np.mean(data[:, i]), "+/-",
        print np.std(data[:, i])
    plt.plot(data[:, 0], data[:, W], 'r-', lw=2)
    plt.plot(data[:, 0], MM(data[:, 0], mm), 'b--', lw=1)
    plt.plot(data[:, 0], data[:, secondMax(data)], 'g-', lw=1)
    plt.xlabel("time (host generations)", fontsize=FS)
    plt.ylabel("number of individuals", fontsize=FS)
    plt.xticks(fontsize=TFS)
    plt.yticks(fontsize=TFS)

    plt.grid(True)

#    maxx = 10 * np.ceil(data[:, 1::].max() / 10.0)
#    fig2 = plt.figure(2, figsize=(10, 6))
#    fig2.canvas.set_window_title("Distribution of sizes of parasite"
#                                 + " No. %s population" % W)
#    hh = plt.hist(data[:, W], color=(lii, lii, lii, 0.95))
#    hmax = 10 * np.ceil(np.max(hh[0]) / 10.0) + 20.0
#    plt.vlines(mm, 0, hmax, colors='b', linestyles='dashed', lw=2)
#    plt.axis([0, maxx, 0, hmax])
#    plt.xlabel("number of individuals")
#    plt.ylabel("frequency")
#    plt.grid(True)
#    print "Done!"
    plt.show()

if __name__ == "__main__":
    main()
