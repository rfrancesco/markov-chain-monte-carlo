#!/bin/bash

# Run simulations with GNU Parallel
# Run with argument dry-run to execute parallel --dry-run

BINARY="xspectrum"
N_JOBS=3				# number of concurrent jobs
Ns="100 125"   #10 20 40 50 80"     
#Nt="$(echo {225..250..25})"
#BETA="0.2"
#mL=10
#BovermL="0.02"
M="0.1 0.08" #1 0.5 0.25 0.20 0.125"
N_MEASURES=$(echo "10^5 + 25000" | bc)			 
LOGFILE="log_{3}"				# log start/finish to file - TODO: do this with parallel directly 
OUT_NAME="out_${BINARY}/{1}"            # output filename string
ARG="./$BINARY {1} {1} {2} $N_MEASURES ::: $Ns :::+ ${M}"

if [ "$1" == 'dry-run' ]; then
	parallel --dry-run $ARG
	exit 0
fi

parallel -j $N_JOBS --bar --joblog log --results  $OUT_NAME $ARG
