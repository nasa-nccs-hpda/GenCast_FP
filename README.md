# GenCast-FP end-to-end workflow

## New Workflow

```bash
singularity build --sandbox gencast-fp-latest docker://nasanccs/gencast-fp:latest
```

/discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest/

```bash
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP -B /discover/nobackup/jacaraba /discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/view/gencast_fp_cli.py -h
```

```bash
Preprocess
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP -B /discover/nobackup/jacaraba /discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/view/gencast_fp_cli.py preprocess --start_date 2024-12-01 --end_date 2024-12-01 --outdir /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast_run
```

```bash
Predict
singularity exec --env PYTHONPATH=/discover/nobackup/jacaraba/development/GenCast_FP --nv -B /discover/nobackup/jacaraba /discover/nobackup/jacaraba/development/GenCast_FP/container/gencast-fp-latest python /discover/nobackup/jacaraba/development/GenCast_FP/gencast_fp/prediction/predict_gencast.py --date "2024-12-01" --input_dir  /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast_run --out_dir  /discover/nobackup/jacaraba/development/GenCast_FP/tests/gencast_prediction
```

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
```bash
cd <dir_name>/GenCast-FP
bash fm_gencast.sh
```
