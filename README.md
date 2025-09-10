# GenCast-FP end-to-end workflow

This workflow is to generate GenCast predictions with GEOS-FP as inputs. Follow the steps below to set up and run.

---

## 1. Clone the Repository on DISCOVER
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

## 3. Excute the workflow
```bash
cd <dir_name>/GenCast-FP
bash fm_gencast.sh
```
