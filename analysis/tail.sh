#!/bin/bash

NMEASURES=$(cat "$1" | awk 'BEGIN { nmeasures = 0 } !/^#/ { nmeasures++ } END { print nmeasures }')
NKEEP="$2"
NDISCARD=$(( $NMEASURES - $NKEEP ))

awk -v N=$NDISCARD -f 'discard.awk' "$1"
