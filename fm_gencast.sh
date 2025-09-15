#!/bin/bash
#SBATCH --job-name=gencast
#SBATCH --output=gencast_%j.out
#SBATCH --error=gencast_%j.err
#SBATCH --time=00:30:00           # walltime (adjust)
#SBATCH --partition=gpu_a100      # GPU partition (check your cluster name)
#SBATCH --constraint=rome
#SBATCH --gres=gpu:1              # request 1 GPU
#SBATCH --cpus-per-task=10        # CPU cores per task (adjust)
#SBATCH --mem=60G                 # memory (adjust as needed)

# Step 0: parse arguments from command line, set up singularity container
START_DATE=""
END_DATE=""
OUTPUT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"  # by default the output dir is current dir
WORKING_DIR=$PWD

while getopts "s:e:o:h" opt; do
    case "$opt" in
        s)
            START_DATE="$OPTARG"
            ;;
        e)
            END_DATE="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        h)
            echo "Usage: $0 [-s <start_date>] [-e <end_date>] [-o <output_directory>]"
            echo "  -s: Start date"
            echo "  -e: End date"
            echo "  -o: Output directory"
            echo "  -h: Show this help message"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 [-s <start_date>] [-e <end_date>] [-o <output_directory>]"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
        *)
            echo "Usage: $0 [-s <start_date>] [-e <end_date>] [-o <output_directory>]"
            exit 1
            ;;
    esac
done

OUTPUT_DIR=$(readlink -f "$OUTPUT_DIR")

echo "Running GenCast prediction."
echo "Start Date: $START_DATE"
echo "End Date: $END_DATE"
echo "Output Directory: $OUTPUT_DIR"

# Validate arguments
if [[ -z "$START_DATE" ]]; then
    echo "Warning: No start date specified"
fi

if [[ -z "$END_DATE" ]]; then
    echo "Warning: No end date specified"
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "No output directory specified, defaulting to PWD."
fi

# Make logs dir, redirect slurm output
mkdir -pv $OUTPUT_DIR/logs
# Set up trap to move files when script exits
mv gencast_${SLURM_JOB_ID}.out "$OUTPUT_DIR/logs/" 2>/dev/null || true
mv gencast_${SLURM_JOB_ID}.err "$OUTPUT_DIR/logs/" 2>/dev/null || true
exec 1>> "$OUTPUT_DIR/logs/gencast_${SLURM_JOB_ID}.out"
exec 2>> "$OUTPUT_DIR/logs/gencast_${SLURM_JOB_ID}.err"

# Make output dir
mkdir -pv $OUTPUT_DIR/outputs

# Load container
module load singularity
container="/discover/nobackup/projects/QEFM/containers/qefm-core-debian-all-aifs-20250609-sandbox"

# Step 1 - Convert GEOS-FP data to ERA5-like format
echo "======================================================================"
echo "Running preprocessing from date $START_DATE to $END_DATE"
echo "======================================================================"
module load anaconda
source activate base
cd $WORKING_DIR/preprocess
python fp2e5.py \
--outdir ${OUTPUT_DIR}/outputs/FP2E \
--start_date $START_DATE \
--end_date $END_DATE

# Step 2 - Generate forecasts using the preprocessed data
cd "$WORKING_DIR/prediction/FMGenCast/graphcast" || exit 1  # Change to pred dir

# Run 1 pred per day in range inside container
echo "======================================================================"
echo "Creating 15-day rollout predictions from date $START_DATE to $END_DATE"
echo "======================================================================"
mkdir -p "${OUTPUT_DIR}/outputs/GenCastRaw"
module load singularity
for ((current=$(date -d "$START_DATE" +%s); current<=$(date -d "$END_DATE" +%s); current+=86400)); do
    current_date=$(date -d "@$current" +%Y-%m-%d)
    echo "======================================================================"
    echo "Predicting 15 day rollout on: $current_date"
    echo "======================================================================"

    singularity exec --nv -B "$WORKING_DIR","$OUTPUT_DIR" "$container" python3 -m fm_gencast.py \
        --date "$current_date" \
        --input_dir "$OUTPUT_DIR/outputs/FP2E" \
        --out_dir "$OUTPUT_DIR/outputs/GenCastRaw"
done

# Step 3 (optional) - Post-process the forecast outputs
# Compute ensemble mean and convert to NetCDF for each forecast lead time

# Run postprocess inside container
echo "======================================================================"
echo "Postprocessing predictions from date $START_DATE to $END_DATE"
echo "======================================================================"
module load anaconda
source activate base
cd $WORKING_DIR/postprocess
# Loop through each day in the date range
for ((current=$(date -d "$START_DATE" +%s); current<=$(date -d "$END_DATE" +%s); current+=86400)); do
    current_date=$(date -d "@$current" +%Y-%m-%d)

    # Extract year, month, and day from current date
    YEAR=$(echo "$current_date" | cut -d'-' -f1)
    MONTH=$(echo "$current_date" | cut -d'-' -f2)
    DAY=$(echo "$current_date" | cut -d'-' -f3)

    # Run postprocess inside container for current date
    echo "======================================================================"
    echo "Postprocessing predictions for date $current_date"
    echo "======================================================================"
    python gencast_cf.py \
        --geos_dir ${OUTPUT_DIR}/outputs/FP2E \
        --pred_dir ${OUTPUT_DIR}/outputs/GenCastRaw \
        --output_dir ${OUTPUT_DIR}/outputs/CF \
        --year "$YEAR" \
        --month "$MONTH" \
        --day "$DAY"
done
