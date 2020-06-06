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
weight_decay=$6
leadtype=$7

echo 'data='$data
echo 'modelPath='$modelPath
echo 'radius='$radius
echo 'T='$T
echo 'weight_decay='$weight_decay

mkdir $modelPath

logPath=$modelPath'/'$data'_'$radius'_'$T'_'$weight_decay
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
    if [ $leadtype == "0" ]
    then
        time python runRandom.py $data0 $modelPath 0 $radius $T 100  $weight_decay 2 > $out0'_100_2_0' 2>&1 & 
        time python runRandom.py $data0 $modelPath 1 $radius $T 100  $weight_decay 5 > $out0'_100_5_0' 2>&1 & 
        time python runRandom.py $data0 $modelPath 2 $radius $T 250  $weight_decay 2 > $out0'_250_2_0' 2>&1 & 
        time python runRandom.py $data0 $modelPath 3 $radius $T 250  $weight_decay 5 > $out0'_250_5_0' 2>&1 & 
    fi
    if [ $leadtype == "1" ]
    then
        time python runRandom.py $data1 $modelPath 0 $radius $T 100  $weight_decay 2 > $out1'_100_2_1' 2>&1 & 
        time python runRandom.py $data1 $modelPath 1 $radius $T 100  $weight_decay 5 > $out1'_100_5_1' 2>&1 & 
        time python runRandom.py $data1 $modelPath 2 $radius $T 250  $weight_decay 2 > $out1'_250_2_1' 2>&1 & 
        time python runRandom.py $data1 $modelPath 3 $radius $T 250  $weight_decay 5 > $out1'_250_5_1' 2>&1 & 
    fi
    if [ $leadtype == "2" ]
    then
        time python runRandom.py $data2 $modelPath 0 $radius $T 100  $weight_decay 2 > $out2'_100_2_2' 2>&1 & 
        time python runRandom.py $data2 $modelPath 1 $radius $T 100  $weight_decay 5 > $out2'_100_5_2' 2>&1 & 
        time python runRandom.py $data2 $modelPath 2 $radius $T 250  $weight_decay 2 > $out2'_250_2_2' 2>&1 & 
        time python runRandom.py $data2 $modelPath 3 $radius $T 250  $weight_decay 5 > $out2'_250_5_2' 2>&1 & 
    fi   
else
    echo $datatype
    if [ $leadtype == "0" ]
    then
        time python runScaffold.py $data0 $modelPath 0 $radius $T 100  $weight_decay 2 > $out0'_100_2_0' 2>&1 & 
        time python runScaffold.py $data0 $modelPath 1 $radius $T 100  $weight_decay 5 > $out0'_100_5_0' 2>&1 & 
        time python runScaffold.py $data0 $modelPath 2 $radius $T 250  $weight_decay 2 > $out0'_250_2_0' 2>&1 & 
        time python runScaffold.py $data0 $modelPath 3 $radius $T 250  $weight_decay 5 > $out0'_250_5_0' 2>&1 & 
    fi
    if [ $leadtype == "1" ]
    then
        time python runScaffold.py $data1 $modelPath 0 $radius $T 100  $weight_decay 2 > $out1'_100_2_1' 2>&1 & 
        time python runScaffold.py $data1 $modelPath 1 $radius $T 100  $weight_decay 5 > $out1'_100_5_1' 2>&1 & 
        time python runScaffold.py $data1 $modelPath 2 $radius $T 250  $weight_decay 2 > $out1'_250_2_1' 2>&1 & 
        time python runScaffold.py $data1 $modelPath 3 $radius $T 250  $weight_decay 5 > $out1'_250_5_1' 2>&1 & 
    fi
    if [ $leadtype == "2" ]
    then
        time python runScaffold.py $data2 $modelPath 0 $radius $T 100  $weight_decay 2 > $out2'_100_2_2' 2>&1 & 
        time python runScaffold.py $data2 $modelPath 1 $radius $T 100  $weight_decay 5 > $out2'_100_5_2' 2>&1 & 
        time python runScaffold.py $data2 $modelPath 2 $radius $T 250  $weight_decay 2 > $out2'_250_2_2' 2>&1 & 
        time python runScaffold.py $data2 $modelPath 3 $radius $T 250  $weight_decay 5 > $out2'_250_5_2' 2>&1 & 
    fi

#    time python runScaffold.py $data0 $modelPath 1 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out0 2>&1 & 
#    time python runScaffold.py $data1 $modelPath 2 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out1 2>&1 & 
#    time python runScaffold.py $data2 $modelPath 3 $radius $T $fingerprint_dim  $weight_decay $learning_rate > $out2 2>&1 & 
fi
wait
