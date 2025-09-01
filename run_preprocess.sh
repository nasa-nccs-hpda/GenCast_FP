#!/bin/bash

#########################################################
# The intermediate files will be stored in $WORKDIR/test/FP2E
# You can change the path by modifying the --outdir parameter in next step
WORKDIR="$1"
#########################################################

#
# Convert GEOS-FP to ERA-5 coordinate
# outdir: output directory, will serve as input directory for GenCast run
# start_date: start date (YYYY-MM-DD)
# end_date: end date (YYYY-MM-DD)

module load anaconda
source activate base
cd ${WORKDIR}/preprocess
python fp2e5.py \
--outdir ${WORKDIR}/test/FP2E \
--start_date 2024-12-01 \
--end_date 2024-12-02
