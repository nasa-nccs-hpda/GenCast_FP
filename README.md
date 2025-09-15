# GenCast-FP end-to-end workflow

This workflow is to generate GenCast predictions with GEOS-FP as inputs. Follow the steps below to set up and run. The workflow currently only works on DISCOVER filesystems.

## New Workflow

### Downloading the Container

```bash
singularity build --sandbox gencast-fp-latest-fix docker://nasanccs/gencast-fp:latest
```

A version of this container is located at (move later to the project space):

```bash
/discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest
```

### All Arguments Here

```bash
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP -B /discover/nobackup/jacaraba /discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/view/gencast_fp_cli.py -h
```

#### Preprocessing

```bash
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP -B /discover/nobackup/jacaraba /discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/view/gencast_fp_cli.py preprocess --start_date 2024-12-01 --end_date 2024-12-01 --outdir /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast_run
```

#### Predict

#### Sandy for Testing Predict

1. Get salloc session

```bash
salloc --partition=gpu_a100 --constraint=rome --ntasks=10 --gres=gpu:1 --mem-per-gpu=100G -t 10:00:00
```

2. Git Clone (you will need this for now since we are making testing, this wont be needed 
later after we embed it inside the container)

```bash
git clone https://github.com/nasa-nccs-hpda/GenCast_FP --branch deployment-operations
```

3. Run prediction for a single day for now (you change it however you think its appropiate, change the paths
to your username or desired locations, the important portion is the PYTHONPATH stuff)

```bash
Predict
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP --nv -B /discover/nobackup/jacaraba /discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/prediction/predict_gencast.py --start_date 2024-12-01 --end_date 2024-12-01 --input_dir  /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast_run --out_dir  /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast_prediction
```
---

## Previous Workflow

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
```bash
cd <dir_name>/GenCast-FP
bash fm_gencast.sh
```
