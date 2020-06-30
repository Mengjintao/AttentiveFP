#!/bin/bash
#SBATCH --job-name=AttentiveFP     # create a short name for your job
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G       	# memory per cpu-core
#SBATCH --time=150:00:00         # total run time limit (HH:MM:SS)

source activate mpnn

data=$1
modelPath=$2
datatype=$3
radius=$4
T=$5
fingerprint_dim=$6
weight_decay=$7
learning_rate=$8

echo 'data='$data
echo 'modelPath='$modelPath
echo 'radius='$radius
echo 'T='$T
echo 'fingerprint_dim='$fingerprint_dim
echo 'weight_decay='$weight_decay
echo 'learning_rate='$learning_rate

mkdir $modelPath

logPath=$modelPath'/'$data'_'$radius'_'$T'_'$fingerprint_dim'_'$weight_decay'_'$learning_rate
echo 'logPath='$logPath

data0='~/jtmeng/SolCuration/org/'$data'/'$data'_org.csv'
data1='~/jtmeng/SolCuration/clean/'$data'_stand.csv'
data2='~/jtmeng/SolCuration/cure1/'$data'_cure.csv'

out0=$logPath'_org.newout'
out1=$logPath'_cln.newout'
out2=$logPath'_cure1.newout'
echo $data0
echo $data1
echo $data2

echo $out0
echo $out1
echo $out2

type="Random"
if [ $datatype == $type ]
then
    echo $datatype
    time python runRandom.py $data0 $modelPath 1 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out0 2>&1 & 
    time python runRandom.py $data1 $modelPath 2 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out1 2>&1 & 
    time python runRandom.py $data2 $modelPath 3 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out2 2>&1 & 
else
    echo $datatype
    time python runScaffold.py $data0 $modelPath 1 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out0 2>&1 & 
    time python runScaffold.py $data1 $modelPath 2 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out1 2>&1 & 
    time python runScaffold.py $data2 $modelPath 3 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out2 2>&1 & 
fi
wait
