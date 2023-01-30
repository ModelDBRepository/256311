#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue March 30 09:03:49 2021.

@author: spiros
"""
import numpy as np
import os
import sys
import brian2

nrun = int(sys.argv[1])
freq = int(sys.argv[2])

brian2.seed(nrun)

# Create the folder with the output spikes
print(f'RUN: {nrun}')
foldername = f'rate{freq}/run_{nrun}'
os.system(f'mkdir -p -v {foldername}')

N = 1000  # number of 'noise' inputs
time_input = 23000 * brian2.ms  # total time of simulation
rate = freq * brian2.Hz  # rate of the Poisson spike generator
P = brian2.PoissonGroup(N, rates=rate)  # Poisson Group
S = brian2.SpikeMonitor(P)  # Record spikes

# Run the simulation
brian2.run(time_input, report='text', report_period=10 * brian2.second)

# Save the data in txt files
fname = 'noise_'
for s in range(len(S.spike_trains())):
    spiketimes = [round(x/brian2.ms, 1) for x in list(S.spike_trains()[s])]
    np.savetxt(f'{foldername}/{fname}{s}.txt',
               spiketimes,
               fmt='%10.1f', newline='\n')
