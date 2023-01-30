#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 20:45:11 2018.

@author: spiros
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import matplotlib
from scipy.stats import sem


def bar_plots(mydict, metric, path_figs, learning):

    my_list = ['Control', 'No_VIPcells', 'No_VIPCR', 'No_VIPCCK', 'No_VIPPVM',
               'No_VIPNVM', 'No_VIPCRtoOLM', 'No_VIPCRtoBC']
    A = mydict
    A_means = []
    A_sems = []
    for case in my_list:
        A_means.append(np.mean(A[case]))
        A_sems.append(scipy.stats.sem(A[case]))

    plt.figure(1, dpi=300)

    y = A_means
    labels = my_list
    N = len(y)
    x = range(N)

    p1, p2, p3, p4, p5, p6, p7, p8 = plt.bar(x, y, yerr=A_sems)

    p1.set_facecolor('blue')
    p2.set_facecolor('red')
    p3.set_facecolor('green')
    p4.set_facecolor('yellow')
    p5.set_facecolor('lightblue')
    p6.set_facecolor('olive')
    p7.set_facecolor('darkmagenta')
    p8.set_facecolor('darkorange')

    plt.xticks(x, labels, rotation='45')
    plt.ylabel(metric, fontsize=16)
    plt.title(metric)

    plt.savefig(path_figs+'/'+learning+'_'+metric +
                '_barplot.pdf', format='pdf', dpi=300)

    plt.cla()
    plt.clf()
    plt.close()

    # Make Boxplots
    A_list = []
    for case in my_list:
        A_list.append(list(mydict[case]))

    plt.figure(1, dpi=300)

    y = A_list
    labels = my_list
    N = len(y)
    x = range(1, N+1)

    # notch shape box plot
    bplot = plt.boxplot(y, notch=True, vert=True, patch_artist=True,
                        labels=labels)  # will be used to label x-ticks

    # fill with colors
    colors = ['blue', 'red', 'green', 'yellow', 'lightblue',
              'olive', 'darkmagenta', 'darkorange']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    for element in ['fliers', 'means', 'medians', 'caps']:
        plt.setp(bplot[element], color='black')

    plt.xticks(x, labels, rotation='45')
    plt.ylabel(metric, fontsize=16)

    plt.savefig(path_figs+'/'+learning+'_'+metric +
                '_boxplot.pdf', format='pdf', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def bar_plots2(mydict, metric, learning, path_figs, baseline):
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    my_list = ['Control', 'No_VIPcells', 'No_VIPCR', 'No_VIPCCK',
               'No_VIPPVM', 'No_VIPNVM', 'No_VIPCRtoBC', 'No_VIPCRtoOLM']
    A = mydict
    A_means = [np.mean(A['Control']), np.mean(A['No_VIPcells']),
               np.mean(A['No_VIPCR']), np.mean(A['No_VIPCCK']),
               np.mean(A['No_VIPPVM']), np.mean(A['No_VIPNVM']),
               np.mean(A['No_VIPCRtoBC']), np.mean(A['No_VIPCRtoOLM'])]
    A_sems = [sem(A['Control']), sem(A['No_VIPcells']), sem(A['No_VIPCR']),
              sem(A['No_VIPCCK']), sem(A['No_VIPPVM']), sem(A['No_VIPNVM']),
              sem(A['No_VIPCRtoBC']), sem(A['No_VIPCRtoOLM'])]

    plt.figure(1, dpi=300)

    y = A_means
    labels = my_list
    N = len(y)
    x = range(N)

    pControl, pVIPcells, pVIPCR, pVIPCCK, pVIPPVM, pVIPNVM, pVIPtoBC, pVIPtoOLM = plt.bar(
        x, y, yerr=A_sems)

    pControl.set_facecolor('blue')
    pVIPcells.set_facecolor('green')
    pVIPCR.set_facecolor('yellow')
    pVIPCCK.set_facecolor('red')
    pVIPPVM.set_facecolor('lightblue')
    pVIPNVM.set_facecolor('lightgreen')
    pVIPtoBC.set_facecolor('yellowgreen')
    pVIPtoOLM.set_facecolor('darkred')

    plt.axhline(y=baseline, linestyle='--', linewidth=2)

    plt.plot()
    plt.xticks(x, labels)
    plt.ylabel(metric, fontsize=16)
    plt.title(learning)
    plt.ylim([0, 0.4])
    plt.savefig(path_figs+learning+'/'+metric +
                '_barplot.pdf', format='pdf', dpi=300)

    plt.cla()
    plt.clf()
    plt.close()

    # Make Boxplots
    A_list = [list(A['Control']), list(A['No_VIPcells']), list(A['No_VIPCR']),
              list(A['No_VIPCCK']), list(A['No_VIPPVM']), list(A['No_VIPNVM']),
              list(A['No_VIPCRtoBC']), list(A['No_VIPCRtoOLM'])]

    plt.figure(1, dpi=300)

    y = A_list
    labels = my_list
    N = len(y)
    x = range(1, N+1)

    # notch shape box plot
    # will be used to label x-ticks
    bplot = plt.boxplot(y, notch=True, vert=True,
                        patch_artist=True, labels=labels)

    # fill with colors
    colors = ['blue', 'green', 'yellow', 'red',
              'lightblue', 'lightgreen', 'yellowgreen', 'darkred']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    for element in ['fliers', 'means', 'medians', 'caps']:
        plt.setp(bplot[element], color='black')

    plt.xticks(x, labels)
    plt.ylabel(metric, fontsize=16)
    plt.title(learning)

    plt.savefig(path_figs+learning+'/'+metric +
                '_boxplot.pdf', format='pdf', dpi=300)

    plt.cla()
    plt.clf()
    plt.close()
