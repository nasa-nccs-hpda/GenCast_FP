#!/bin/bash

#########################################################
# The intermediate files will be stored in $WORKDIR/test/CF
# You can change the path by modifying the --output_dir parameter
WORKDIR="$1"
#########################################################

#
# Convert GEOS-FP to ERA-5 coordinate
# outdir: output directory, will serve as input directory for GenCast run
# start_date: start date (YYYY-MM-DD)
# end_date: end date (YYYY-MM-DD)

module load anaconda
source activate base
cd ${WORKDIR}/postprocess
python gencast_cf.py \
--geos_dir ${WORKDIR}/test/FP2E \
--pred_dir ${WORKDIR}/test/GenCastRaw \
--output_dir ${WORKDIR}/test/CF \
--year "2024" \
--month "12" \
--day "01"
