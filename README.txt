
# README file for Shuman et al, 2019 Breakdown of spatial coding and interneuron synchronization in epileptic mice. Nature Neuroscience

# for more information, refer to the comments inside the scripts or contact me in: chavlis [DOT] spiros [AT] gmail [DOT] com

# Scripts' author: S. Chavlis, PhD, I. Pandi, M.Sc.


#################### INPUT CREATION #################################################################################
# First you have to create the Inputs, go in make_inputs_linear_track directory
cd make_inputs_linear_track

In a command line execute

python make_grid_like_inputs.py <run_number>
python sp_make_place_inputs.py <run_number>
python glim_v2_prelearning.py <run_number> <desynch_level> <jitter_source>

# <run_number> is a specific run from one edge of the track to the other. To replicate the figures one needs 10 runs
# <desynch_level>: ms of desynchronization
# <jitter_source>: which input to randomize. Valid options: EC or CA3. Use EC to replicate the paper figures.

#e.g., python glim_v2_prelearning.py 1 20 EC

# Then enter background_noise directory
cd ../background_noise

# create the background noise by executing 

python poisson_input.py <total_number_of_runs> <poisson_rate>  # e.g., poisson_input.py 1 5 --> creates run1 poisson random noise with lambda 5 Hz


#################### MAIN SIMULATIONS #######################################################################
# return to main directory
cd ../

# Compile all mechanisms (mod files)
nrnivmodl mechanisms/

# Run the simulation
.x86_64/special -nogui -c nruns=<run_number> -c ntrials=<virtual_mouse_id> -c desynch=<desynch_level_in_ms> -c n_neuron=<deletion_type> -c factor=<reduction_factor> Network_CA1.hoc

#e.g., .x86_64/special -nogui -c nruns=1 -c ntrials=1 -c n_neuron=0 -c desynch=0 -c factor=1 Network_CA1.hoc

# to replicate the results of the paper you need 10 runs/trial and 10 trials and all possible deletions (see below)

# Valid deletions: 
# Control: All connections and cells, default  -- option: 0
# SOMred:  SOMs are removed by a specific percentage  -- option: 1
# PVred:   PVs are removed by a specific percentage  -- option: 2
# Desynch: Desynchronization of EC/CA3 inputs by a specific amount  -- option: 3
# ALL:     SOMred, PVred and Desynch simultaneously  -- option: 4
# SOMdel:  All PVs are removed  -- option: 5
# PVdel:   All PVs are removed  -- option: 6


# Output of the simulation is saved into Simulation_Results/

#################### ANALYSIS OF LOCOMOTION DATA BEFORE PROCEEDING #####################################################

# First, one needs to extract the spiketimes for neurons in order to analyze them
# Go to AnalysisRawData directory
cd AnalysisRawData

# Exctract spike times

python spiketimes_analysis.py <neuron_type> <deletion_type> <number_of_trial> <number_of_run>

# Valid deletions_types: 
# Control
# SOMred
# PVred
# Desynch
# ALL
# SOMdel
# PVdel

# Valid <neuron_type> values:

# _pvsoma_   : Pyramidal cells
# _aacell_   : Axoaxonic cells
# _bcell_    : Basket cells
# _bscell_   : Bistratified cells
# _olmcell_  : OLM cells
# _vipcck_   : VIP/CCK cells
# _vipcr_    : VIP/CR PVM cells
# _vipcrnvm_ : VIP/CR NVM cells

# After the analysis for all trials, runs and deletions execute:

python all_path_all_spiketimes.py <deletion_type> # e.g., python all_path_all_spiketimes.py 0

# This will create the subfolder final_results/metrics_permutations
# where the spiketimes and the path for all cases is stored (for better handling)

# Permutations for all cells to find spatial information and stability null distributions
python permutations_analysis.py <virtual_mouse> <pyramidalID> <deletion_type>

# Data save for using in GraphPad Prism and basic plotting
python analysis_path.py <virtual_mouse> <deletion_type>

python all_trials_paper_all.py

python all_trials_per_animal.py

# for more information, refer to the comments inside the scripts or contact me in: chavlis [DOT] spiros [AT] gmail [DOT] com







