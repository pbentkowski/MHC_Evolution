#!/bin/bash

# Script for launching the MHC Evolution Model on the PCSS Eagle cluster that is
# using Slurm queueing system. The script will submit as many jobs to the
# Slurm queuing system as many lines are in the ParamParam.csv file. Each simulation
# will write its results to a separate directory.

if [ $# -ne 2 ]
  then
    echo "Give the name of the simulation's executable, please."
    echo "And the number j to enumerate directories. The numbering will start from j+1."
    exit 1
fi

IFS=$'\n'
j=$2 # numbering of directories will start from j+1
for line in $(< ParamParam.csv);
  do
    j=$(($j+1))
    mkdir MHC.$j
    cp $1 MHC.$j/
    cd MHC.$j
    IFS=$' '
    # the lines below will stitch a custom sbatch Slurm file and run it
    echo "srun ./$1 $line" > theLine.txt
    cat ~/launchers/sbatchTemplate_begin.sh theLine.txt ~/launchers/sbatchTemplate_end.sh > runTheMHCjob.sl
    rm theLine.txt
    sbatch --job-name=$j.mhc runTheMHCjob.sl &
    echo -e "Run No. $j launched! Params are set to: $line"
    IFS=$'\n'
    cd ..
  done
  wait
w=$(($((`squeue -u $USER | wc -l`))-1))
echo -e "We have $w jobs launched!"
echo -e "Ride'em computing nodes!"

