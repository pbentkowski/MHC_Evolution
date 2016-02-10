#! /bin/bash


# Launches a defined number of program instances (see $procc) placing each in its own
# directory. The script feeds the simulations with arguments taken from ParamParam.csv file.
# Total number of program runs is equal to the number of lines in file ParamParam.csv,
# but the number of programs run simultaneously in never larger than the number specified
# by the $procc variable.

if [ $# -eq 0 ]
  then
    echo "Give the name of the model's executable, please."
    exit 1
fi

if [ -x $1 ]
 then
    procc=4  # number of cores/processors used for calculations
    j=0 # numbering of directories will start from j+1
    execdir=$PWD  # model's exec file has to be in the same directory as this script
    IFS=$'\n'
    for line in $(< ParamParam.csv);
      do
        while true
          do
            modelNumb=$(($((`ps -x | grep $1 | wc -l`))-3))
      if (( $modelNumb < $procc )) ; then
        j=$(($j+1))
        mkdir MHC.$j
        cd MHC.$j
        IFS=$' '
        nohup $execdir/$1 $line &
        echo -e "Run No. $j launched! Params are set to: $line"
        IFS=$'\n'
        cd ..
        break
      else
        sleep 15m
      fi
        done
      done;
    echo -e "All modells launched!"
  else
    echo -e "Sorry, but $1 is not an executable file :-( "
fi
