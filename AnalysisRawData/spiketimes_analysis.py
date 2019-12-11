import sys,pickle,time
from peakdet import peakdet
import numpy as np

#Loads a list with path points for each time-series point

time1=time.time()

def extract_spiketimes(neuron_type,condition,ntrial,nrun):
    '''
    Extract Spike Times from NEUORN's voltage trace
    Nneurons: number of neurons to be analyzed
    condition: Specific Lesion
    ntrial: mouse ID
    nrun: number of trial/run
    '''
    dt = 0.1
    filepath = '../Simulation_Results/'+condition+'/Trial_'+str(ntrial)+'/Run_'+str(nrun)
    
    # Peakdet parameters
    delta = 1
    
    if neuron_type=='_pvsoma_':
        Nneurons=130
        thres = 0
    elif neuron_type=='_bcell_':
        Nneurons=8
        thres = 0
    elif neuron_type=='_vipcck_' or neuron_type=='_vipcrnvm_':
        Nneurons=1
        thres = 0
    elif neuron_type=='_vipcr_':
        Nneurons=4
        thres = 0        
    else:
        Nneurons=2
        thres = 0
        
    print Nneurons
    spiketimes_all=[]
    for n_neuron in xrange(Nneurons):

        
        filename = filepath+'/Trial_'+ntrial+'_Run_'+nrun+neuron_type+str(n_neuron)+'.dat'
    
        path = np.loadtxt('../make_inputs_linear_track/runs_produced_by_python_ec_rand_stops/run_'
                          +nrun+'/path.txt','int',delimiter=' ')
        tot_time = np.sum(np.bincount(path[:,0])) # in ms
        
        data = np.loadtxt(filename)
        # remove the first 400ms
        data = data[int(400/dt):int(tot_time/dt)+1]

        maxtab, mintab = peakdet(data, delta, thres)
        if maxtab.size!=0:    
            spiketimes = [round(float(i*dt),1) for i in maxtab[:,0]]
        else:
            spiketimes=[]
        spiketimes.insert(0,str(n_neuron))
        
        spiketimes_all.append(spiketimes)
    
    #Saves the list in a pickle file at the specified Run directory
    filewrite = filepath + '/spiketimes'+ neuron_type +'.pkl'
    with open(filewrite, 'wb') as handle:
        pickle.dump(spiketimes_all, handle, protocol=pickle.HIGHEST_PROTOCOL)
            

ntype    = sys.argv[1]
cond     = sys.argv[2]
trial    = sys.argv[3]
run      = sys.argv[4]

results = extract_spiketimes(ntype, cond, trial, run)

print "\nEverything was ok for pyramidals. Case: "+cond+" "+str(trial)+" "+str(run)
time2=time.time()
duration = round(time2-time1,3)
print "\n The analysis run for "+str(duration)+" seconds"
