#/bin/bash

# you can set this to any directory you want
WORKDIR=$PWD

# Step 1 - Convert GEOS-FP data to ERA5-like format
bash $WORKDIR/run_preprocess.sh $WORKDIR && echo "Preprocessing done"
if [ $? -ne 0 ]; then
    echo "Preprocessing failed"
    exit 1
fi

# Step 2 - Generate forecasts using the preprocessed data
# Create logs directory if it doesn't exist
mkdir -p $WORKDIR/logs
sbatch --wait $WORKDIR/run_prediction.sh $WORKDIR && echo "Forecasting job submitted"

# Step 3 (optional) - Post-process the forecast outputs
# Compute ensemble mean and convert to NetCDF for each forecast lead time
bash $WORKDIR/run_postprocess.sh $WORKDIR && echo "Postprocessing done"
if [ $? -ne 0 ]; then
    echo "Postprocessing failed"
    exit 1
fi
