#!/bin/bash

if [ -z "$1" ]; then
  echo "[info] No detector file specified. Using default: ATOMX.detector"
  DETECTOR="ATOMX.detector"
else
  DETECTOR="$1"
fi

DETNAME="${DETECTOR%.detector}"
OUTPUTNAME="simulation_${DETNAME}"
echo "$OUTPUTNAME"

if [ -z "$2" ]; then
  RANDOMSEED=20250517
else
  RANDOMSEED="$2"
fi

npsimulation --detector "$DETECTOR" \
             --event-generator r_34Ar_alpha_p.reaction\
             --random-seed "$RANDOMSEED"\
             --output "$OUTPUTNAME"\
             -B geant4_run.mac
