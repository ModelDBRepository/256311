#--------------------------------------------------- Batch 1


read -p "Enter 'EC' or 'CA3' jitter --> " value

counter=0
for nruns in $(seq 1 50); do
for desynch in $(seq 0 5 50);do
	export 	LPARAMS="${nruns} ${desynch} ${value}"
	echo $LPARAMS
	qsub -v "LPARAMS=$LPARAMS" submit_spiros_final_inputs.sh
	counter=$[counter+1]
done
done
