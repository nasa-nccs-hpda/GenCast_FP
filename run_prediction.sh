#!/bin/bash
#SBATCH --job-name=gencast_pred
#SBATCH --output=logs/gencast_pred_%j.out
#SBATCH --error=logs/gencast_pred_%j.err
#SBATCH --time=00:30:00           # walltime (adjust)
#SBATCH --partition=gpu_a100      # GPU partition (check your cluster name)
#SBATCH --constraint=rome
#SBATCH --gres=gpu:1              # request 1 GPU
#SBATCH --cpus-per-task=10        # CPU cores per task (adjust)
#SBATCH --mem=60G                 # memory (adjust as needed)

# Load necessary modules
module load singularity

# Container path
#container="/discover/nobackup/projects/QEFM/containers/qefm-core-gencast-20250511-sandbox"
container="/discover/nobackup/projects/QEFM/containers/qefm-core-debian-all-aifs-20250609-sandbox"

# Working directory
WORKDIR="/discover/nobackup/jli30/GenCast_FP"

cd "$WORKDIR/prediction/FMGenCast/graphcast" || exit 1

# Run inside container
singularity exec --nv -B "$WORKDIR" "$container" python3 -m fm_gencast.py


# # Load necessary modules

# module load singularity
# #container="/discover/nobackup/projects/QEFM/containers/qefm-core-gencast-20250511-sandbox"
# container="/discover/nobackup/projects/QEFM/containers/qefm-core-debian-all-aifs-20250609-sandbox"

# # Define working directory and input files
# WORKDIR="/discover/nobackup/jli30/GenCast_FP"

# cd "$WORKDIR/prediction/FMGenCast/graphcast" || exit 1

# singularity exec -B "$WORKDIR" "$container" python3 -m fm_gencast.py
