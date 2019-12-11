#--------------------------------------------------- Batch 1
counter=0
for n_trials in $(seq 1 10); do
for n_runs in $(seq 1 10); do
	export 	LPARAMS="-c n_trials=${n_trials} -c n_runs=${n_runs} -c n_neuron=$1 -c desynch=$2 -c factor=$3"
	echo $LPARAMS
	qsub -v "LPARAMS=$LPARAMS" submit_spiros_peyman_plc3.sh
	counter=$[counter+1]
done
done
