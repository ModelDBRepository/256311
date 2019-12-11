# declare a name for this job to be sample_job

#Maximum walltime for this job
#$ -l h_rt=9999:00:00
#Maximum cpu time for this job
#$ -l h_cpu=9999:00:00

# Specify the shell to use when running the job script
#$ -S /bin/sh

# Directory to perform the job
#$ -cwd

# Name of the Job
#$ -N CA1Epil

#$ -o output_files
#$ -e output_files

# run the program
./x86_64/special -nogui $LPARAMS Network_CA1.hoc
