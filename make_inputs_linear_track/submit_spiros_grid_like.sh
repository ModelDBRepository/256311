# declare a name for this job to be sample_job

# request 1 node
#PBS -l nodes=1:ppn=1
#PBS -l cput=9999:00:00
#PBS -l walltime=240:00:00

#$ -S /bin/sh
#$ -cwd
#$ -N GridLike
#$ -j n

#$ -o log_files_cluster/make_GridInputs.$JOB_ID.$TASK_ID.out
#$ -e log_files_cluster/make_GridInputs.$JOB_ID.$TASK_ID.err

echo "Starting job: $SGE_TASK_ID"

# run the program
python make_grid_like_inputs.py $SGE_TASK_ID

echo "Done with job: $SGE_TASK_ID"
