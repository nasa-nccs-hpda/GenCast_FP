#!/bin/bash

## NOTE; Modify the paths below to match your environment
WORKDIR='/discover/nobackup/jli30/GenCast_FP'
#########################################################

#
# step 1
# Convert GEOS-FP to ERA-5 coordinate
#

module load anaconda
source activate base
cd ${WORKDIR}/preprocess
python fp2e5.py \
--outdir ${WORKDIR}/output_test \
--start_date 2024-12-01 \
--end_date 2024-12-02
