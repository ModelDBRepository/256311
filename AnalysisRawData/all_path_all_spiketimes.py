#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 15:18:16 2018.

@author: spiros
"""
import sys
import numpy as np
import os
import pickle


def spiketimes_path_load(arg1):
    case = arg1
    fsave = 'final_results/metrics_permutations/'
    os.system('mkdir -p '+fsave)

    my_list = [case]
    nruns = 10
    trials = 10
    Npyramidals = 130
    for ntrial in range(1, trials+1):
        print
        print(f'TRIAL: {ntrial}')
        print
        spiketimes_all = {}
        path_all = []
        for nrun in range(1, nruns+1):
            print(f'RUN: {nrun}')

            # Load PATH and SPIKETIMES of Pyramidal npyr
            # Load of path -- different for each run
            pathd = '../make_inputs_linear_track/runs_produced_by_python_ec_rand_stops/run_' + \
                str(nrun)
            path = np.loadtxt(pathd+'/path.txt', 'int', delimiter=' ')
            path_all.append(path)

            Aspiketimesdict = {}
            for npyr in range(Npyramidals):
                Aspiketimes = {}
                for case in my_list:
                    cond = '../Simulation_Results'
                    fileload = cond+'/'+case+'/Trial_' + \
                        str(ntrial)+'/Run_'+str(nrun)+'/spiketimes_pvsoma_.pkl'

                    if os.path.isfile(fileload):
                        with open(fileload, 'rb') as f:
                            spiketimes_load = pickle.load(f)
                        spiketimes = spiketimes_load[npyr]
                        # remove first entry -- aka pyramidal number
                        spiketimes = spiketimes[1:]
                        spiketimes = [
                            x for x in spiketimes if x < path.shape[0]]
                    else:
                        spiketimes = []
                        print("File does not exist.")
                        continue

                    Aspiketimes[case] = spiketimes
                Aspiketimesdict['Pyramidal'+str(npyr)] = Aspiketimes
            spiketimes_all['Run'+str(nrun)] = Aspiketimesdict

        with open(fsave+'path_all_trial_'+str(ntrial)+'_'+case+'.pkl', 'wb') as handle:
            pickle.dump(path_all, handle, protocol=pickle.HIGHEST_PROTOCOL)

        with open(fsave+'spiketimes_all_trial_'+str(ntrial)+'_'+case+'.pkl', 'wb') as handle:
            pickle.dump(spiketimes_all, handle,
                        protocol=pickle.HIGHEST_PROTOCOL)


arg1 = sys.argv[1]
spiketimes_path_load(arg1)
