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
sbatch $WORKDIR/run_prediction.sh $WORKDIR && echo "Forecasting job submitted"

