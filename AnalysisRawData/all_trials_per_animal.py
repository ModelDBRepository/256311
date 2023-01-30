#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:01:48 2017.

@author: spiros
"""

import pandas as pd
import pickle
import os
import scipy.stats

import matplotlib.pyplot as plt
from place_cell_metrics import sparsity_index2
from place_cell_metrics import selectivity_index
from place_cell_metrics import peak_frequency
from place_cell_metrics import field_size, spatial_coherence
import numpy as np

# matplotlib.use('agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from visualization_ import make_dicts, bar_plots

# Save in MATLAB format
from scipy.io import savemat

npath_x, npath_y = 200, 1
Nbins = 100
trials = 10
Npyramidals = 130
Nperms = 500

what_to_do = 3


if what_to_do == 0:
    my_list = ['Control', 'ALL0.50_20', 'Desynch20', 'SOMred0.50', 'PVred0.50']
elif what_to_do == 1:
    my_list = ['Control']
    for i in range(5, 36, 5):
        my_list += ['Desynch'+str(i), 'ALL0.50_'+str(i)]
elif what_to_do == 2:
    my_list = ['Control']
    for i in ['0.75', '0.50', '0.25']:
        my_list += ['SOMred'+i, 'PVred'+i]
    my_list += ['SOMdel', 'PVdel']
elif what_to_do == 3:
    my_list = ['Control']
    for i in range(5, 36, 5):
        my_list += ['Desynch'+str(i), 'ALL0.50_'+str(i)]
    for i in ['0.75', '0.50', '0.25']:
        my_list += ['SOMred'+i, 'PVred'+i]
    my_list += ['SOMdel', 'PVdel']

spec = 'final_results/'
path_figs = spec+'/figures2/'
path_figs2 = spec+'/figures2_metrics/'
file_load = spec+'/metrics2/'
file_load_perms = spec+'/metrics_permutations/'
path_data = spec+'/data_final3/'

trials = [str(i) for i in range(1, trials+1)]

os.system('mkdir -p ' + path_figs2)
os.system('mkdir -p '+path_data)

numbers_all = {}
information_all = {}
information_plc = {}
stability_index_plc = {}
stability_index_all = {}
stability_all_plc = {}
stability_all_all = {}
sparsity_all = {}
selectivity_all = {}
peak_freq_all = {}
mean_freq_in_all = {}
mean_freq_out_all = {}
fieldsize_all = {}
fieldsize_plc = {}
coherence_all = {}
coherence_plc = {}

SpatialMapsALL = {}
SpatialMapsPLC = {}

RateMapMat = {}

for case in my_list:
    rate_map_mat = np.zeros((len(trials), Nbins, Npyramidals))
    for trial in trials:
        with open(file_load+'pickled_sn_'+case+'_'+trial+'.pkl', 'rb') as f:
            loaded_data = pickle.load(f)
        for npyr in range(Npyramidals):
            rate_map = loaded_data['maps'].squeeze()
            rate_map = rate_map[npyr, :]
            rate_map_mat[int(trial)-1, :, int(npyr)] = rate_map

    RateMapMat[case] = rate_map_mat

for ntrial in trials:
    print("TRIAL:", ntrial)
    rateMaps = {}
    time_bin = {}
    information1 = {}
    information2 = {}
    stability_index = {}
    sparsity = {}
    selectivity = {}
    peak_freq = {}
    mean_freq_in = {}
    mean_freq_out = {}
    fieldsize1 = {}
    fieldsize2 = {}
    numbers_plc = {}
    number_of_peaks = {}
    reward_zone = {}
    stability1 = {}
    stability2 = {}
    stability1_all = {}
    stability2_all = {}
    coherence1 = {}
    coherence2 = {}

    SpatialMapsALL['Mouse'+str(ntrial)] = {}
    SpatialMapsPLC['Mouse'+str(ntrial)] = {}

    for case in my_list:
        with open(file_load+'/pickled_sn_'+case+'_'+ntrial+'.pkl', 'rb') as f:
            loaded_data = pickle.load(f)

        place_cells = np.loadtxt('../Simulation_Results/Control/Trial_' +
                                 str(ntrial)+'/Run_1/input_plcs.txt', delimiter=',')
        place_cells = [int(x) for x in list(place_cells[:, 0])]
        rateMaps[case] = loaded_data['maps']
        time_bin[case] = loaded_data['time_in_bin']

        inforALL1 = []
        stabALL1 = []
        stabALL11 = []
        inforALL2 = []
        stabALL2 = []
        stabALL22 = []
        sparsALL = []
        selecALL = []
        peaksALL = []
        sizesALL1 = []
        sizesALL2 = []
        average1ALL = []
        average2ALL = []
        numbersALL = 0
        place_cells_idx = []
        place_cells_max = []
        coher1 = []
        coher2 = []

        for npyr in range(Npyramidals):

            with open(file_load_perms+'perms_pickled_info_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                infor = pickle.load(f)
            with open(file_load_perms+'perms_pickled_stab_even_odd_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                stab1 = pickle.load(f)
            with open(file_load_perms+'perms_pickled_stab_first_second_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                stab2 = pickle.load(f)
            with open(file_load_perms+'perms_pickled_stab_all_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                stab_all = pickle.load(f)
            stab_all = [np.nanmean(i) for i in stab_all]

            pval_infor = sum(infor[0] >= infor[1:])/float(Nperms)
            pval_stab1 = sum(stab1[0] >= stab1[1:])/float(Nperms)
            pval_stab2 = sum(stab2[0] >= stab2[1:])/float(Nperms)
            pval_staball = sum(stab_all[0] >= stab_all[1:])/float(Nperms)

            rate_map = rateMaps[case][npyr, :, :]
            t_bin = time_bin[case][npyr, :, :]

            maxpeak = np.max(rate_map)
            meanrate = np.mean(rate_map)
            sizetest1 = field_size(
                rate_map, relfreq=0.2*maxpeak, track_length=Nbins)[0]
            infield = field_size(rate_map, relfreq=0.2 *
                                 maxpeak, track_length=Nbins)[1]
            outfield = field_size(rate_map, relfreq=0.2 *
                                  maxpeak, track_length=Nbins)[2]
            sp_coher = spatial_coherence(rate_map.squeeze(), window=3)
            clevel = 0.99
            if (maxpeak >= 1.0) and (pval_infor >= clevel) and (pval_staball >= clevel) and (5./(npath_x/Nbins) <= sizetest1 <= 25./(npath_x/Nbins)):

                numbersALL += 1
                place_cells_idx.append(npyr)
                place_cells_max.append(np.argmax(rate_map))

                inforALL2.append(infor[0])
                stabALL2.append(
                    (np.math.atanh(stab1[0])+np.math.atanh(stab2[0]))/2.0)
                stabALL22.append(np.math.atanh(stab_all[0]))
                sparsALL.append(sparsity_index2(rate_map))
                selecALL.append(selectivity_index(rate_map))
                peaksALL.append(peak_frequency(rate_map))
                sizesALL2.append(sizetest1)
                average1ALL.append(field_size(
                    rate_map, relfreq=0.2*np.max(rate_map), track_length=Nbins)[1])
                average2ALL.append(field_size(
                    rate_map, relfreq=0.2*np.max(rate_map), track_length=Nbins)[2])
                coher2.append(sp_coher)

            inforALL1.append(infor[0])
            stabALL1.append(
                (np.math.atanh(stab1[0])+np.math.atanh(stab2[0]))/2.0)
            stabALL11.append(np.math.atanh(stab_all[0]))
            sizesALL1.append(sizetest1)
            coher1.append(sp_coher)

        rate_maps_all = rateMaps[case][:, :, :]
        Ncells = rate_maps_all.shape[0]

        idx = np.argmax(rate_maps_all.squeeze(), axis=1)
        new_idx = np.lexsort((range(Ncells), idx))
        rtMaps = rate_maps_all[new_idx, :, :].squeeze()

        Max = np.max(rtMaps, axis=1).reshape(-1, 1)
        for i in range(Max.shape[0]):
            if Max[i, 0] == 0:
                Max[i, 0] = 1e-12

        rtMaps = rtMaps / np.repeat(Max, Nbins, axis=1)

        fig, axes = plt.subplots(nrows=1, ncols=2)
        SpatialMapsALL['Mouse'+str(ntrial)][case] = rtMaps
        im0 = axes[0].imshow(rtMaps.squeeze(), cmap="jet", aspect='equal')
        divider = make_axes_locatable(axes[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im0, cax=cax)
        ax = plt.gca()

        np.savetxt(file_load+'place_cell_idx_trial_'+str(ntrial)+'_'+case+'.txt',
                   place_cells_idx, fmt='%i')

        rate_maps_plc = rateMaps[case][place_cells_idx, :, :]
        Ncells = rate_maps_plc.shape[0]

        idx = np.argmax(rate_maps_plc.squeeze(), axis=1)
        new_idx = np.lexsort((range(Ncells), idx))
        rtMaps = rate_maps_plc[new_idx, :, :].squeeze()

        Max = np.max(rtMaps, axis=1).reshape(-1, 1)
        for i in range(Max.shape[0]):
            if Max[i, 0] == 0:
                Max[i, 0] = 1e-12
        rtMaps = rtMaps / np.repeat(Max, Nbins, axis=1)

        SpatialMapsPLC['Mouse'+str(ntrial)][case] = rtMaps
        im1 = axes[1].imshow(rtMaps.squeeze(), cmap="jet", aspect='equal')
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(axes[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im1, cax=cax)
        plt.tight_layout()
        plt.savefig(path_figs+'/'+case+'_PlaceCells_heatmap_' +
                    str(ntrial)+'.pdf', format='pdf', dpi=300)
        plt.cla()
        plt.clf()
        plt.close()

        print(numbersALL/130., case)
        numbers_plc[case] = numbersALL
        information1[case] = inforALL1
        stability1[case] = stabALL1
        stability1_all[case] = stabALL11
        information2[case] = inforALL2
        stability2[case] = stabALL2
        stability2_all[case] = stabALL22
        sparsity[case] = sparsALL
        selectivity[case] = selecALL
        peak_freq[case] = peaksALL
        mean_freq_in[case] = average1ALL
        mean_freq_out[case] = average2ALL
        fieldsize1[case] = sizesALL1
        fieldsize2[case] = sizesALL2
        coherence1[case] = coher1
        coherence2[case] = coher2

        if case in numbers_all.keys():
            numbers_all[case].append(numbers_plc[case]/130.)
        else:
            numbers_all[case] = [numbers_plc[case]/130.]

        information_all = make_dicts(information_all, information1, case)
        information_plc = make_dicts(information_plc, information2, case)
        stability_index_all = make_dicts(stability_index_all, stability1, case)
        stability_index_plc = make_dicts(stability_index_plc, stability2, case)
        stability_all_all = make_dicts(stability_all_all, stability1_all, case)
        stability_all_plc = make_dicts(stability_all_plc, stability2_all, case)
        sparsity_all = make_dicts(sparsity_all, sparsity, case)
        selectivity_all = make_dicts(selectivity_all, selectivity, case)
        peak_freq_all = make_dicts(peak_freq_all, peak_freq, case)
        mean_freq_in_all = make_dicts(mean_freq_in_all, mean_freq_in, case)
        mean_freq_out_all = make_dicts(mean_freq_out_all, mean_freq_out, case)
        fieldsize_all = make_dicts(fieldsize_all, fieldsize1, case)
        fieldsize_plc = make_dicts(fieldsize_plc, fieldsize2, case)
        coherence_all = make_dicts(coherence_all, coherence1, case)
        coherence_plc = make_dicts(coherence_plc, coherence2, case)

mydict_all = {}
mydict_all['numbers_place_cells'] = numbers_all
mydict_all['information_all'] = information_all
mydict_all['stabilityIndex_all'] = stability_index_all
mydict_all['stabilityALL_all'] = stability_all_all
mydict_all['information_plc'] = information_plc
mydict_all['stabilityIndex_plc'] = stability_index_plc
mydict_all['stabilityALL_plc'] = stability_all_plc
mydict_all['sparsity'] = sparsity_all
mydict_all['selectivity'] = selectivity_all
mydict_all['fieldsize_all'] = fieldsize_all
mydict_all['fieldsize_plc'] = fieldsize_plc
mydict_all['peak_freq'] = peak_freq_all
mydict_all['numbers_plc'] = numbers_all
mydict_all['mean_in'] = mean_freq_in_all
mydict_all['mean_out'] = mean_freq_out_all
mydict_all['coherence_all'] = coherence_all
mydict_all['coherence_plc'] = coherence_plc


mydict_all1 = {}
mydict_all1['Description'] = '10 animals, 10 trials/animal, 130 pyramidal cells/animal'
mydict_all1['number_of_place_cells'] = numbers_all
mydict_all1['information_of_all_cells'] = information_all
mydict_all1['information_of_place_cells'] = information_plc
mydict_all1['stability_of_all_cells'] = stability_index_all
mydict_all1['stability_of_place_cells'] = stability_index_plc
mydict_all1['spatial_maps_of_all_cells'] = SpatialMapsALL
mydict_all1['spatial_maps_of_place_cells'] = SpatialMapsPLC

savemat(path_data+'Data_spiros.mat', mydict_all1, oned_as='column')

# placecells = numbers_all.values()
placecells = []
for cased in my_list:
    placecells.append(numbers_all[cased])

fig = plt.figure(1, figsize=(5, 5))
y = placecells
y = [i for i in list(np.mean(placecells, axis=1))]
ye = [i for i in list(scipy.stats.sem(placecells, axis=1))]
labels = my_list
N = len(y)
x = range(N)

colors = ['blue', 'red', 'green', 'yellow',
          'lightblue', 'olive', 'darkmagenta', 'darkorange']
plt.bar(x, y, color=colors[:len(my_list)], yerr=ye)

plt.xticks(x, labels, rotation='90')
plt.ylabel('Number of Place Cells', fontsize=16)
plt.savefig(path_figs2+'/'+'numberOfPlaceCells_barplot.eps',
            format='eps', dpi=300)
plt.savefig(path_figs2+'/'+'numberOfPlaceCells_barplot.png',
            format='png', dpi=300)
plt.cla()
plt.clf()
plt.close()


placecells = []
for cased in my_list:
    placecells.append(numbers_all[cased])

plt.figure(1, dpi=150)

y = placecells
labels = my_list
N = len(y)
x = range(1, N+1)

# notch shape box plot
bplot = plt.boxplot(y, notch=False, vert=True, patch_artist=True,
                    labels=labels)  # will be used to label x-ticks

# fill with colors
colors = ['blue', 'red', 'green', 'yellow',
          'lightblue', 'olive', 'darkmagenta', 'darkorange']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

for element in ['fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot[element], color='black')

plt.xticks(x, labels, rotation='45')
plt.ylabel('number of place cells', fontsize=16)
plt.savefig(path_figs2+'/'+'numberOfPlaceCells_boxplot.eps',
            format='eps', dpi=300)
plt.savefig(path_figs2+'/'+'numberOfPlaceCells_boxplot.png',
            format='png', dpi=300)
plt.cla()
plt.clf()
plt.close()


my_list2 = ['information_all', 'information_plc', 'sparsity',
            'selectivity', 'peak_freq', 'fieldsize_all', 'fieldsize_plc', 'stabilityIndex_all',
            'stabilityIndex_plc', 'mean_in', 'mean_out', 'stabilityALL_all',
            'stabilityALL_plc', 'coherence_all', 'coherence_plc']

if what_to_do == 0:
    fnam = '_figure_pval_'
elif what_to_do == 1:
    fnam = '_Desynch_pval_'
elif what_to_do == 2:
    fnam = '_INs_pval_'
elif what_to_do == 3:
    fnam = '_ALL_pval_'

for metric in my_list2:
    bar_plots(mydict_all[metric], metric, path_figs2, my_list)


# USE this for PRISM-GraphPad plotting
my_list2 = ['information_all', 'information_plc', 'sparsity',
            'selectivity', 'peak_freq', 'fieldsize_all', 'fieldsize_plc', 'stabilityIndex_all',
            'stabilityIndex_plc', 'mean_in', 'mean_out', 'stabilityALL_all',
            'stabilityALL_plc', 'coherence_all', 'coherence_plc']

L = len(my_list)
for metric in my_list2:
    B = np.zeros((len(trials), L*3))
    cnt = 0
    for case in my_list:
        path_data = spec+'/data_final3/'
        A = np.array(mydict_all[metric][case])
        B[:, 3*cnt:3*(cnt+1)] = A
        cnt += 1

    df = pd.DataFrame(B)
    df.to_csv(path_data+metric+fnam+str(clevel)+'.csv', sep=' ', index=False)
