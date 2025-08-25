#!/bin/bash

# This script runs the prediction using the GenCast model and GEOS-FP input data.

# Load necessary modules

module load singularity
container="/discover/nobackup/projects/QEFM/containers/qefm-core-gencast-20250507-sandbox"

# Define working directory and input files
WORKDIR="/discover/nobackup/jli30/GenCast_FP/prediction/FMGenCast"
cd "$WORKDIR" || exit 1

singularity exec "$container" python3 -m fm_gencast.py
