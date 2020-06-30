#!/bin/bash
batch=(64 32 16)
hidsize=(9 11 13)

data=$1
modelPath=$2
datatype=$3

for radius in 2 4 6
do
	for T in 1 3 5
        do
		for fingerprint_dim in 100 250
                do
                        for weight_decay in 2 4 6 
                        do
                        	for learning_rate in 2 5 
                        	do
					sbatch -p szsc szbatch.sh $data $modelPath $datatype $radius $T $fingerprint_dim $weight_decay $learning_rate 	
				#	$sh szbatch.sh $data $modelPath $radius $T $fingerprint_dim $weight_decay $learning_rate
				done
			done
         	done
         done
done
