#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue March 30 09:05:21 2021.

@author: spiros
"""
import os
import sys
import numpy as np

nrun = int(sys.argv[1])
desync_factor = sys.argv[2]
jitter = sys.argv[3]

np.random.seed(nrun)

npath_x, npath_y = 200, 1

time_delay = 400
time_delay_ca3 = 17

n_place_fields = 41

dend_ec = 8
dend_ca3 = 6

if jitter == 'EC':
    name = 'EC'
elif jitter == 'CA3':
    name = 'CA3'
else:
    sys.exit('Wrong jittering.')

folderpath = f'Inputs_linear_rand_stops_noisy_jitter{name}_{desync_factor}/'

if not os.path.exists(folderpath):
    os.system(f'mkdir -p {folderpath}')

rndAll = []

vec_shuffled_ec = []
vec_shuffled_ca3 = []

print(f'Factor {desync_factor} RUN {nrun} ...\n')

p = f'run_{nrun}/'
source_ec = f'runs_produced_by_python_ec_rand_stops/{p}'
source_ca3 = f'runs_produced_by_python_ca3_rand_stops/{p}'
dest = folderpath + p

if not os.path.exists(dest):
    os.makedirs(dest)

listdirs = []
counter_ec = 0
counter_ca3 = 0

for plf in range(1, n_place_fields+1):
    listdirs.append(f'place_field_{plf}')

for mydir in listdirs:

    L_ec = len(os.listdir(source_ec+mydir))
    vec_ec = list(np.random.permutation(range(dend_ec)))
    vec_shuffled_ec.append(vec_ec)

    for i in range(L_ec):

        with open(f'{source_ec}{mydir}/s{vec_ec[i]}.txt', 'r') as f:
            lines = f.read().splitlines()

        lines = [int(x)+time_delay for x in lines]
        lines = sorted(lines)
        maximum = lines[-1]
        lines_noisy = []

        for iline in lines:
            if np.random.rand() <= 0.10:
                lines_noisy.append(int(maximum*np.random.rand()))
            else:
                lines_noisy.append(iline)

        lines_noisy = list(set(lines_noisy))
        if jitter == 'EC':
            jitter = list(np.round(int(desync_factor) *
                                   np.random.randn(len(lines_noisy)), 1))
            lines_noisy = list(set([x+y for x, y in zip(lines_noisy, jitter)]))

        lines_noisy = sorted(lines_noisy)
        lines_noisy = [x for x in lines_noisy if x > 0]

        with open(f'{dest}g{counter_ec}_EC.txt', 'w') as f:
            for line in lines_noisy:
                nline = str(int(line))
                f.write(nline + '\n')

        counter_ec += 1

print('ok with ec')

for mydir in listdirs:

    L_ca3 = len(os.listdir(source_ca3+mydir))
    vec_ca3 = list(np.random.permutation(range(dend_ca3)))
    vec_shuffled_ca3.append(vec_ca3)

    for i in range(0, L_ca3):

        with open(f'{source_ca3}{mydir}/c{vec_ca3[i]}.txt', 'r') as f:
            lines = f.read().splitlines()

        lines = [int(x) + time_delay+time_delay_ca3 for x in lines]
        lines = sorted(lines)
        maximum = lines[-1]
        lines_noisy = []

        for iline in lines:
            if np.random.rand() <= 0.05:
                lines_noisy.append(int(maximum*np.random.rand()))
            else:
                lines_noisy.append(iline)

        lines_noisy = list(set(lines_noisy))
        if jitter == 'CA3':
            jitter = list(np.round(int(desync_factor) *
                                   np.random.randn(len(lines_noisy)), 1))
            lines_noisy = list(set([x+y for x, y in zip(lines_noisy, jitter)]))

        lines_noisy = sorted(lines_noisy)
        lines_noisy = [x for x in lines_noisy if x > 0]

        with open(f'{dest}g{counter_ca3}_CA3.txt', 'w') as f:
            for line in lines_noisy:
                nline = str(int(line))
                f.write(nline + '\n')

        counter_ca3 += 1

print('ok with ca3')
