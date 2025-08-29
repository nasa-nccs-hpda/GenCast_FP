#!/bin/bash

# This script runs the prediction using the GenCast model and GEOS-FP input data.

# Load necessary modules

module load singularity
#container="/discover/nobackup/projects/QEFM/containers/qefm-core-gencast-20250511-sandbox"
container="/discover/nobackup/projects/QEFM/containers/qefm-core-debian-all-aifs-20250609-sandbox

# Define working directory and input files
WORKDIR="/discover/nobackup/jli30/GenCast_FP"

cd "$WORKDIR/prediction/FMGenCast/graphcast" || exit 1

singularity exec -B "$WORKDIR" "$container" python3 -m fm_gencast.py
