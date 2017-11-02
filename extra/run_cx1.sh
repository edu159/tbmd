#!/bin/sh 
#PBS -l walltime=20:00:00 
#PBS -l select=1:ncpus=12:mem=20gb

module load gsl
module load intel-suite
path=$(pwd)
input="input"
output="output"
config="config"

exec_path="$path/tbmd -i $input -o $output -c $config"
echo $exec_path
#eval $exec_path
