# GenCast-FP end-to-end workflow

This workflow is to generate GenCast predictions with GEOS-FP as inputs. Follow the steps below to set up and run.

---

## 1. Clone the Repository
```bash
mkdir $WORK_DIR
cd $WORK_DIR
git clone https://github.com/nasa-nccs-hpda/GenCast_FP.git
cd $WORK_DIR/GenCast-FP
```

## 2. Copy checkpoints and ancillary dataset
```bash
cd $WORK_DIR/GenCast-FP/prediction
cp /discover/nobackup/jli30/GenCast_FP/prediction/checkpoint.tar.gz .
tar -xzvf checkpoint.tar.gz
cd ..
```

## 3. Excute the workflow
```bash
cd $WORK_DIR/GenCast-FP
bash fm_gencast.sh
```
