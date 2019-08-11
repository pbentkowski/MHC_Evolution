#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "Give a number how many lines from top of the file to keep?"
    exit 1
fi

echo "$(head -n $1 HostsGeneDivers.csv)" > HostsGeneDivers.csv
echo "$(head -n $1 NoMutationInPathoList.csv)" > NoMutationInPathoList.csv
echo "$(head -n $1 NumberOfMhcAfterMating.csv)" > NumberOfMhcAfterMating.csv
echo "$(head -n $1 NumberOfMhcBeforeMating.csv)" > NumberOfMhcBeforeMating.csv
echo "$(head -n $1 NumberOfMhcInFather.csv)" > NumberOfMhcInFather.csv
echo "$(head -n $1 NumberOfMhcInMother.csv)" > NumberOfMhcInMother.csv
echo "$(head -n $1 PresentedPathogenNumbers.csv)" > PresentedPathogenNumbers.csv
