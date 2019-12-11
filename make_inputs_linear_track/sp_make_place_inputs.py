import os, sys
import numpy as np

nrun    = int(sys.argv[1])

maindir = 'runs_produced_by_python_ca3_rand_stops' ### ###prelearning 1 (io makes the ca3 input) #neg pos the same input
source = 'runs_produced_by_python_ec_rand_stops'  ###prelearning 1 (io makes the ec input) #neg pos the same input

if not os.path.exists(maindir):
    os.system('mkdir -p '+maindir)

np.random.seed(nrun)

nplace_field = 41
step         = 5
nEC          = 8
nCA3         = 6
npath_x      = 200
npath_y      = 1

half_size    = 6

# Place like inputs parameters
mu_in  = 0.16
mu_out = 0.016
sigma  = 0.002

print 'RUN:...' + str(nrun)
folder = maindir + '/run_' + str(nrun)
os.system('mkdir -p '+ str(folder))

############################################################################
######## when ec inputs(.txt) store in matlab without time map ##########
pathd = source+'/run_'+str(nrun)
path = np.loadtxt(pathd+'/path.txt','int',delimiter=' ')

path=path[:,0]

L = len(path)

time_array = np.bincount(path)[1:]

csum = np.cumsum(time_array)
csum = np.insert(csum,0,0)

xarray = range(0,npath_x+1,step) # plf begins from 1 to nplace_field+1
yarray = [1]


for plf in range(1, nplace_field+1): #for plf 1 --> restriction if
    print "Place Field... ", str(plf)

    folder2 = folder + '/place_field_' + str(plf)
    os.system('mkdir -p '+ str(folder2))

    # Load ALL spiketimes/octal
    spikemap_sall=[]
    for dend in range(0, nEC):
        with open(source+'/run_'+str(nrun)+'/place_field_'+str(plf)+'/s'+str(dend)+'.txt','r', 0) as f:
            A = f.read().splitlines()
            spikemap_sall+=[int(x) for x in A]

    vector = sorted(spikemap_sall)
    peak = xarray[plf-1]  # plf in range(1,nplace_field+1): begins from 1..

###################################################################
####################### CHECK THE PLFS  #############################
#####################################################################

    #########################
    initial = peak-half_size
    final   = peak+half_size

    # Check boundary conditions
    if initial < 0:
        initial=0

    if final > npath_x:
        final=npath_x
    ########################     
    
    inplf_ca3_input  = []
    outplf_ca3_input = []

    for spiketime in vector:
        if spiketime>=csum[initial] and spiketime<=csum[final]:                
            inplf_ca3_input.append(spiketime)
        else:
            outplf_ca3_input.append(spiketime)


    for i in range(0,nCA3):        

        z_in   = (mu_in + np.random.randn()*sigma)
        z_out  = (mu_out + np.random.randn()*sigma)

        shuf_vec_in  = []
        shuf_vec_out = []
        
        inplf_ca3_input = list(np.random.permutation(inplf_ca3_input))
        outplf_ca3_input = list(np.random.permutation(outplf_ca3_input))
        
        count_in  = int(len(inplf_ca3_input) * z_in)
        count_out = int(len(outplf_ca3_input) * z_out)
        
        # Boundary condition
        if count_in==0:
            count_in=-1
        if count_out==0:
            count_out=-1
        #####################

        shuf_vec_in  = inplf_ca3_input[-count_in:]
        shuf_vec_out = outplf_ca3_input[-count_out:]

        ca3_input = shuf_vec_in+shuf_vec_out

        ca3_input = sorted(set(ca3_input))
        
        ca3_input = np.array(ca3_input).reshape(-1,1)
        filename = folder2+'/c'+str(i)+'.txt'
        np.savetxt(filename, ca3_input, fmt='%d', delimiter=' ')
