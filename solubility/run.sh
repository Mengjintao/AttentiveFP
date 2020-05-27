#!/bin/bash
#SBATCH --job-name=AttentiveFP     # create a short name for your job
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G       	# memory per cpu-core
#SBATCH --time=150:00:00         # total run time limit (HH:MM:SS)

source activate mpnn

data=$1
modelPath=$2
mkdir $modelPath

data0='~/jtmeng/SolCuration/org/'$data'/'$data'_org.csv'
data1='~/jtmeng/SolCuration/clean/'$data'_stand.csv'
data2='~/jtmeng/SolCuration/cure1/'$data'_cure.csv'

out0=$modelPath'/'$data'_org.newout'
out1=$modelPath'/'$data'_cln.newout'
out2=$modelPath'/'$data'_cure1.newout'
echo $data0
echo $data1
echo $data2

echo $out0
echo $out1
echo $out2

time python run.py $data0 $modelPath 1 > $out0 2>&1 & 
time python run.py $data1 $modelPath 2 > $out1 2>&1 & 
time python run.py $data2 $modelPath 3 > $out2 2>&1 & 
wait
 
