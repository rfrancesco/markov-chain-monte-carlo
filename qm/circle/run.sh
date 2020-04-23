#!/bin/bash

# Run simulations with GNU Parallel
# Run with argument dry-run to execute parallel --dry-run

ALGORITHM=tailor
BINARY="x${ALGORITHM}"
N_JOBS=4				# number of concurrent jobs
NETAS="5"
NS="200 300"
N_MEASURES=$(echo "10^5" | bc)			 
N_SKIP=1
LOGFILE=log				# log start/finish to file - TODO: do this with parallel directly 
OUT_NAME="out_${ALGORITHM}/{1}/{2}"            # output filename string
ARG="./$BINARY {1} {2} $N_MEASURES $N_SKIP ::: $NETAS ::: ${NS}"

if [ "$1" == 'dry-run' ]; then
	parallel --dry-run $ARG
	exit 0
fi

echo "$BINARY: $(date): Started job: Ns = [${NS}], Neta = [$NETAS]" >> $LOGFILE
parallel -j $N_JOBS --bar --results $OUT_NAME $ARG
echo "$BINARY: $(date): Finished job: Ns = [${NS}], Neta = [$NETAS]" >> $LOGFILE
