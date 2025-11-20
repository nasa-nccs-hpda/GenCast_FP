# GenCast-FP end-to-end workflow

[![DOI](https://zenodo.org/badge/1042752062.svg)](https://doi.org/10.5281/zenodo.17662736)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/nasa-nccs-hpda/GenCast_FP/README.md)
[![Build Docker](https://github.com/nasa-nccs-hpda/GenCast_FP/actions/workflows/dockerhub.yml/badge.svg?event=release)](https://github.com/nasa-nccs-hpda/GenCast_FP/actions/workflows/dockerhub.yml)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/nasa-nccs-hpda/GenCast_FP)
![Docker Image Version](https://img.shields.io/docker/v/nasanccs/gencast-fp?label=Docker)
![License](https://img.shields.io/github/license/nasa-nccs-hpda/GenCast_FP)

This workflow is to generate GenCast predictions with GEOS-FP as inputs.
Follow the steps below to set up and run. The workflow currently only works on DISCOVER filesystems.

## Quickstart

The following command runs preprocessing, prediction, and postprocessing for a given date
range using the Discover A100 systems. You will need access to a single GPU to run this workflow.
Note that the following command can be run from any Discover login node.

For a single day (end_date defaults to the same day):

```bash
sbatch --partition=gpu_a100 --constraint=rome --ntasks=10 --gres=gpu:1 \
    --mem-per-gpu=100G -t 10:00:00 -J gencast-fp \
    --wrap="module load singularity; singularity exec --nv \
    -B $NOBACKUP,/css,/gpfsm/dmd/css,/nfs3m,/gpfsm \
    /discover/nobackup/projects/QEFM/containers/gencast-fp-latest \
    gencast-fp run --start_date 2025-11-20:00 \
    --output_dir /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast-run"
```

if you want to run for multiple past days:

```bash
sbatch --partition=gpu_a100 --constraint=rome --ntasks=10 --gres=gpu:1 \
    --mem-per-gpu=100G -t 10:00:00 -J gencast-fp \
    --wrap="module load singularity; singularity exec --nv \
    -B $NOBACKUP,/css,/gpfsm/dmd/css,/nfs3m,/gpfsm \
    /discover/nobackup/projects/QEFM/containers/gencast-fp-latest \
    gencast-fp run --start_date 2025-11-10:00 --end_date 2025-11-15:00 \
    --output_dir /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast-run"
```

Example slurm file submission script:

```bash
#!/bin/bash
#SBATCH --partition=gpu_a100
#SBATCH --constraint=rome
#SBATCH --ntasks=10
#SBATCH --gres=gpu:1
#SBATCH --mem-per-gpu=100G
#SBATCH --time=10:00:00
#SBATCH --job-name=gencast-fp
#SBATCH --output=gencast-fp_%j.out
#SBATCH --error=gencast-fp_%j.err


# Load modules
source /usr/share/modules/init/bash
module purge
module load singularity

# Run the container
singularity exec --nv \
    -B $NOBACKUP,/css,/gpfsm/dmd/css,/nfs3m,/gpfsm \
    /discover/nobackup/projects/QEFM/containers/gencast-fp-latest \
    gencast-fp run \
    --start_date 2025-11-20:00 \
    --output_dir /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast-run
```

## Dependencies

Additional details and flexibility of the commands are listed below.

### Downloading the Container

If you would like to download the container yourself, you will need to run the following
command. The latest version has the most up to date changes, while specific releases are
attached to a given version from the repository.

#### Latest Release

```bash
singularity build --sandbox gencast-fp-latest docker://nasanccs/gencast-fp:latest
```

#### Specific Version

```bash
singularity build --sandbox gencast-fp-0.1.0 docker://nasanccs/gencast-fp:0.1.0
```

A version of this container is located at:

```bash
/discover/nobackup/projects/QEFM/containers/gencast-fp-latest
```

## Pipeline Details

In addition, individual steps of the pipeline can be run using the container and CLI. Some examples with arguments
are listed below. The pipeline has 3 steps: preprocess, predict, and postprocess.

### Preprocessing

```bash
singularity exec --nv -B $NOBACKUP,/css,/gpfsm/dmd/css,/nfs3m,/gpfsm /discover/nobackup/projects/QEFM/containers/gencast-fp-latest_jordan gencast-fp preprocess -h
WARNING: underlay of /usr/bin/nvidia-smi required more than 50 (225) bind mounts
usage: gencast_fp_cli.py preprocess [-h] --start_date START_DATE --end_date END_DATE [--output_dir OUTPUT_DIR] [--expid EXPID] [--res_value RES_VALUE] [--nsteps NSTEPS]

options:
  -h, --help            show this help message and exit
  --start_date START_DATE
                        Start date to process (YYYY-MM-DD:HH)
  --end_date END_DATE   End date to process (YYYY-MM-DD:HH)
  --output_dir OUTPUT_DIR
                        Output directory for preprocessed files
  --expid EXPID         Experiment ID for the output files
  --res_value RES_VALUE
                        Resoluton (default 1.0 resolution)
  --nsteps NSTEPS       Number of steps for rollout (default 30, 15 days)
```

### Prediction

```bash
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP --nv -B $NOBACKUP,/css,/gpfsm/dmd/css,/nfs3m,/gpfsm /discover/nobackup/projects/QEFM/containers/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/view/gencast_fp_cli.py predict -h
usage: gencast_fp_cli.py predict [-h] --start_date START_DATE --end_date END_DATE --input_dir INPUT_DIR --output_dir OUTPUT_DIR [--ckpt CKPT] [--nsteps NSTEPS] [--res RES] [--ensemble ENSEMBLE]
                                 [--container_meta CONTAINER_META]

options:
  -h, --help            show this help message and exit
  --start_date START_DATE
                        YYYY-MM-DD
  --end_date END_DATE   YYYY-MM-DD
  --input_dir INPUT_DIR, -i INPUT_DIR
                        Preprocessed input directory
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Where to write predictions
  --ckpt CKPT           Path to GenCast .npz checkpoint (overrides container default)
  --nsteps NSTEPS
  --res RES
  --ensemble ENSEMBLE
  --container_meta CONTAINER_META
                        Where to load default ckpt/configs if --ckpt not passed
```

### Postprocessing

```bash
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP --nv -B $NOBACKUP,/css,/gpfsm/dmd/css,/nfs3m,/gpfsm /discover/nobackup/projects/QEFM/containers/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/view/gencast_fp_cli.py postprocess -h
usage: gencast_fp_cli.py postprocess [-h] --start_date START_DATE --end_date END_DATE --input_dir INPUT_DIR --predictions_dir PREDICTIONS_DIR [--output_dir OUTPUT_DIR] [--no_ens_mean]

options:
  -h, --help            show this help message and exit
  --start_date START_DATE
                        Start date (YYYY-MM-DD)
  --end_date END_DATE   End date (YYYY-MM-DD)
  --input_dir INPUT_DIR
                        Directory with GEOS inputs (for initial conditions)
  --predictions_dir PREDICTIONS_DIR
                        Directory with GenCast predictions
  --output_dir OUTPUT_DIR
                        Directory for CF-compliant NetCDF outputs
  --no_ens_mean         Disable ensemble mean (keep all ensemble members)
```
