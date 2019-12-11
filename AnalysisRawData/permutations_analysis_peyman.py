#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 21:04:08 2018

@author: spiros
"""
import sys, os, time
from permutations_shuffle import permutations_shuffle

def permutations(ntrial, npyr, case):
    import pickle
    print "Analyze: " + case
    print
    
    folder1='final_results/metrics_permutations/'
    os.system('mkdir -p '+folder1)
    
    fileload = 'final_results/metrics_permutations/'
    
    with open(fileload+'path_all_trial_'+str(ntrial)+'_'+case+'.pkl', 'rb') as f:
        path_all=pickle.load(f)    
    
    with open(fileload+'spiketimes_all_trial_'+str(ntrial)+'_'+case+'.pkl', 'rb') as f:
        spiketimes_all=pickle.load(f)
    
    print "Pyramidal: " + str(npyr)
    mydict = permutations_shuffle(case, npyr, path_all, spiketimes_all,npath_x=200, 
                                  npath_y=1,Nperms=501, Nbins=100, nruns=10)
        
    infop              = mydict['information']     
    stability_even_odd = mydict['stabEvenOdd']
    stability_fir_sec  = mydict['stabFirstSecond']
    stability_all      = mydict['stabAll']
    rateMaps_all       = mydict['rateMaps']
    
    
              
    with open(folder1+'perms_pickled_info_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+str(ntrial)+'.pkl', 'wb') as handle:
        pickle.dump(infop, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(folder1+'perms_pickled_stab_even_odd_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+str(ntrial)+'.pkl', 'wb') as handle:
        pickle.dump(stability_even_odd, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(folder1+'perms_pickled_stab_first_second_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+str(ntrial)+'.pkl', 'wb') as handle:
        pickle.dump(stability_fir_sec, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(folder1+'perms_pickled_stab_all_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+str(ntrial)+'.pkl', 'wb') as handle:
        pickle.dump(stability_all, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(folder1+'perms_pickled_rateMap_all_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+str(ntrial)+'.pkl', 'wb') as handle:
        pickle.dump(rateMaps_all, handle, protocol=pickle.HIGHEST_PROTOCOL)        

tic      = time.time()
ntrial   = sys.argv[1]
npyr     = sys.argv[2] 
case     = sys.argv[3]
results  = permutations(ntrial, npyr, case )
toc      = time.time()

print "\nTotal time: "+str(round(toc-tic,3))+" seconds"

