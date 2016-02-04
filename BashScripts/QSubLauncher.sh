#! /bin/bash

# Script for launching the MHC Evolution Model on the PL-Grid clusters that are
# using PBS queuing system. File ParamParam.csv contains in each subsequent line
# parameters fed to the model run. The script will submit as many jobs to the
# PBS queuing system as many lines in ParamParam.csv there are. Each simulation
# will write its results to a separate directory.
# See (in Polish): https://docs.plgrid.pl/pages/viewpage.action?pageId=4260614

IFS=$'\n'
j=0 # numbering of directories will start from j+1
jobName='MHC_evo'
for line in $(< ParamParam.csv);
  do
    j=$(($j+1))
    mkdir MHC.$j
    cp mhcevolution MHC.$j/
    cd MHC.$j
    IFS=$' '
    echo "cd $PWD && $PWD/mhcevolution $line" | qsub -N $jobName \
		-l walltime=70:23:23 -l mem=2gb -A plgpbentkowsk2015a
    echo -e "Run No. $j launched! Params are set to: $line"
    IFS=$'\n'
    cd ..
  done;
w=$(($((`qstat -u $USER | wc -l`))-5))
echo -e "We have $w modells launched!"
