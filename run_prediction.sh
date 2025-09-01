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

# Create logs directory if it doesn't exist
mkdir -p logs

# Load necessary modules
module load singularity

# Container path
container="/discover/nobackup/projects/QEFM/containers/qefm-core-debian-all-aifs-20250609-sandbox"

# Working directory setting in "fm_gencast.sh"
WORKDIR=$1
echo "Working directory: $WORKDIR"

cd "$WORKDIR/prediction/FMGenCast/graphcast" || exit 1

# Run inside container
singularity exec --nv -B "$WORKDIR" "$container" python3 -m fm_gencast.py \
--date "2024-12-01" \
--input_dir "$WORKDIR/test/FP2E" \
--out_dir "$WORKDIR/test/GenCastRaw" \
