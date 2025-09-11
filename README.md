# GenCast-FP end-to-end workflow

This workflow is to generate GenCast predictions with GEOS-FP as inputs. Follow the steps below to set up and run. The workflow currently only works on DISCOVER filesystems.

---

## 1. Clone the Repository on DISCOVER
**Note: ensure that there is ample space (>100 gb) to run this code. DISCOVER nobackup is ideal.**
```bash
mkdir <dir_name>
cd <dir_name>
git clone https://github.com/nasa-nccs-hpda/GenCast_FP.git
cd <dir_name>/GenCast-FP
```

## 2. Copy checkpoints and ancillary dataset
```bash
cd <dir_name>/GenCast-FP/prediction
cp /discover/nobackup/jli30/GenCast_FP/prediction/checkpoint.tar.gz .
tar -xzvf checkpoint.tar.gz
cd ..
```

## 3. Execute the workflow
The wrapper script needs to be run using sbatch to utilize GPU resources. The script makes one 15-day rollout prediction per day in the range: [start_date, end_date], inclusive. The arguments to the script are as follows: 

* -s: start date of predictions, in YYYY-MM-DD format (dashes are required)
* -e: end date of predictions, in YYYY-MM-DD format (dashes are required)
* -o: output directory, defaults to current working directory if not specified. Files will be placed under a subdir of this directory, called outputs/. Three directories will be housed in outputs:
  * FP2E: GEOS-FP data converted to ERA5, output of preprocessing step.
  * GenCastRaw: Raw prediction NetCDF files, representing one 15-day rollout per file.
  * CF: Postprocessed 15-day rollouts, 1 rollout per 12-hour slice (30 outputs in total). Postprocessed data is separated into subdirectories based on year, month, and day of the start of the rollout (so a rollout on 2024-12-01 will have 30 files in a subfolder with a name like: Y2024/M12/D01). Postprocessed data contains ensemble mean.

Below is an example, using a start date of 12/01/2024 and an end date of 12/02/2024 (thus 2 15-day rollouts will be generated). 
```bash
cd <dir_name>/GenCast-FP
sbatch --wait fm_gencast.sh -s 2024-12-01 -e 2024-12-02 -o $PWD
```
