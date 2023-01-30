#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 09:08:09 2018.

@author: spiros
"""
import os
import sys
import numpy as np
from gridfield import gridfield
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap

# Adopted from https://github.com/BIDS/colormap/blob/master/parula.py
cm_data = np.loadtxt('parula_like_colormap.txt')
parula_map = LinearSegmentedColormap.from_list('parula', cm_data)

# Choose coordinates +1 for the value
# e.g., if you want 100-->myx = 101
# Apart from the case that one dimension is one


##############################################################################
# P A T H   C O N S T R U C T I O N  #########################################
##############################################################################
my_run = int(sys.argv[1])

visualize = False  # if True, uncomment line 37, comment 36

myx = 201  # Choose +1 the points you want (apart from 1)
myy = 1

# Place field coordinations; all to all combinations
x_array = range(0, myx, 5)
y_array = [1]

print(f'Simulating RUN... {my_run}')
print
maindir = 'runs_produced_by_python_ec_rand_stops/'

dirname = f'{maindir}run_{my_run}'
os.system(f'mkdir -p {dirname}')

np.random.seed(my_run)

nx = 1
ny = 1

p0 = [nx, ny]
k = 1

zold = 0

totlength = 0
if myy == 1:
    npoints = (myx-1)*myy
else:
    npoints = (myx-1)*(myy-1)

nstops = np.random.randint(0, 10, 1).item()
print(f'Number of Random Stops: {nstops}')
rand_stops = sorted(list(set(np.random.randint(1, myx-2, nstops))))
rand_stops_minus = [i-1 for i in rand_stops]
rand_stops_plus = [i+1 for i in rand_stops]

path = []
for i in range(1, npoints+1):
    if (ny > myy):
        break

    if (i in rand_stops) or (i in rand_stops_minus) or (i in rand_stops_plus):
        sigma = 50
        mu = 200
        print('Random Stop')
    else:
        sigma = 30
        mu = 50
    # random time at each point
    # mu = 50
    # random time at each point
    z = int(mu + np.random.randn(1)*sigma)

    while (z < 1):
        z = int(mu + np.random.randn(1)*sigma)

    # for same time at each point
    # z=mu;

    time_at_the_same_point = z
    path += [p0]*int(time_at_the_same_point)

    nx = nx + k
    p0 = [nx, ny]

    if (nx > myx-1) and (ny % 2) != 0:
        ny = ny+1
        nx = myx-1
        p0 = [nx, ny]
        k = -1

    if (nx < 1) and (ny % 2) == 0:
        ny = ny+1
        nx = 1
        p0 = [nx, ny]
        k = 1

    zold += time_at_the_same_point

# save the path
path = np.array(path)
filename = f'{dirname}/path.txt'
np.savetxt(filename, path, fmt='%.0f', delimiter=' ')

print('Done with the path')

##############################################################################
#  G R I D    L I K E    I N P U T S  ########################################
##############################################################################
ndend = 8  # Number of dendrites
theta_freq = 8  # in Hz
theta_phase = 0
my_field = 0

for xxx in x_array:
    for yyy in y_array:

        my_field += 1
        folder2 = f'{dirname}/place_field_{my_field}'
        os.system(f'mkdir -p {folder2}')

        # d is the x,y point of the grid field of dend ni
        d = np.zeros((ndend, myx, myy))
        dd = np.zeros((myx, myy))

        angle = 0.0
        lambda_var = 3.0
        for ni in range(ndend):
            lambda_var += 0.5
            angle += 0.4
            for x in range(myx):
                # d is the point x,y of grid field of dend ni
                for y in range(myy):
                    d[ni, x, y] = gridfield(angle, lambda_var, xxx, yyy, x, y)

        for ni in range(ndend):
            dd += d[ni, :, :]

        # run visualize only for few dendrites, is RAM consuming
        if visualize:
            cmap = parula_map

            fig, axes = plt.subplots(nrows=4, ncols=2, dpi=150)
            ni = 0
            for ax in axes.flat:
                im = ax.imshow(d[ni, :, :].T, origin='lower',
                               cmap=cmap, aspect='50', vmin=0, vmax=1)
                ax.tick_params(axis='y', which='both',
                               right='off', left='off', labelleft='off')
                ax.set_xlabel('Position (cm)')
                ax.set_xticks(range(0, 201, 50))
                ax.set_xticklabels([str(x) for x in range(0, 201, 50)])
                ax.set_title('Grid-like Input '+str(ni+1))
                ni += 1

            cbar = fig.colorbar(im, ax=axes.ravel().tolist())
            cbar.ax.set_ylabel('Normalised firing rate [Hz]')
            cbar.set_ticks([x/10.0 for x in range(0, 11, 2)])
            cbar.set_ticklabels([str(x/10.0) for x in range(0, 11, 2)])

            plt.tight_layout()
            plt.savefig('TheoreticalGridCell.eps', format='eps', dpi=1200)
            plt.savefig('TheoreticalGridCell.png', format='png', dpi=1200)

            fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150)
            im = plt.imshow(dd.T/np.max(dd), origin='lower',
                            cmap=cmap, aspect='50', vmin=0, vmax=1)
            ax.tick_params(axis='y', which='both', right='off',
                           left='off', labelleft='off')
            ax.set_xlabel('Position (cm)')
            ax.set_xticks(range(0, 201, 50))
            ax.set_xticklabels([str(x) for x in range(0, 201, 25)])
            ax.set_title('Place-like Cell')
            cbar = plt.colorbar(im)
            cbar.ax.set_ylabel('Normalised firing rate [Hz]')
            plt.savefig('TheoreticalPlaceCell.eps', format='eps', dpi=1200)
            plt.savefig('TheoreticalPlaceCell.png', format='png', dpi=1200)

        for ni in range(ndend):
            spikes = []
            for i in range(len(path)):  # i is time
                current_loc = path[i, :]

                probability = d[ni, current_loc[0]-1, current_loc[1]-1]
                probability *= (np.sin(2.0*np.pi*theta_freq *
                                       i/1000.0 + theta_phase)+1.0)/2.0

                r_ = np.random.rand(1)
                if (probability > 0.7) and (r_ < probability / 2.0):

                    # spikes is a vector with the locations/timespots
                    # where there is a spike
                    spikes.append(i)

            spikes = np.array(spikes).reshape(-1, 1)
            filename2 = f'{folder2}/s{ni}.txt'
            np.savetxt(filename2, spikes, fmt='%.0d', delimiter=' ')

        print(f'Done with Grid field {my_field}')
