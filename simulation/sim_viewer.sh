#!/bin/bash

if [ -z "$1" ]; then
  echo "[info] No detector file specified. Using default: ATOMX.detector"
  DETECTOR="ATOMX.detector"
else
  DETECTOR="$1"
fi

DETNAME="${DETECTOR%.detector}"
OUTPUTNAME="viewer_${DETNAME}"
echo "$OUTPUTNAME"

npsimulation --detector "$DETECTOR" \
             --event-generator r_34Ar_alpha_p.reaction \
             --random-seed 15 \
             --output "$OUTPUTNAME"\
             -M geant4_viewer.mac
