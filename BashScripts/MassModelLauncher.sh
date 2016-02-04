#! /bin/bash


# Launches a defined number of program instances (see $procc) placing each in its own
# directory. The script gives to the simulations arguments taken from ParamParam.csv file.
# Total number of program runs is equal to the number of lines in file ParamParam.csv,
# but the number of programs run simultaneously in never larger than the number specified
# by the $procc variable.


IFS=$'\n'
procc=6  # number of cores/processors used for calculations
j=0 # numbering of directories will start from j+1
for line in $(< ParamParam.csv);
  do
      while true
      do
	       modNo=$(($((`ps -x | grep mhcevolution | wc -l`))-1))
	if (( $modNo < $procc )) ; then
	  j=$(($j+1))
	  mkdir MHC.$j
	  cd MHC.$j
	  echo -e "Run No. $j launched! Params are set to: $line"
	  IFS=$' '
	  /home/piotr/Documents/Projects/mhcEvolution/dist/Release/GNU-Linux-x86/mhcevolution $line &
	  IFS=$'\n'
	  cd ..
	  break
	else
	  sleep 15m
	fi
      done
  done;
echo -e "All modells launched!"
