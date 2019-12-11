# declare a name for this job to be sample_job

#Maximum walltime for this job
#$ -l h_rt=9999:00:00
#Maximum cpu time for this job
#$ -l h_cpu=9999:00:00
#$ -l h_vmem=8G

# Remove a specific computational machine
##$ -l h=!compute-0-18
##$ -l h=!compute-0-13

# Specify the shell to use when running the job script
#$ -S /bin/sh

# Directory to perform the job
#$ -cwd

# Name of the Job
#$ -N finalInputs

#$ -o log_files_cluster/
#$ -e log_files_cluster/

# run the program
python glim_shuf_new_noisy.py $LPARAMS
