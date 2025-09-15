FROM python:3.10-slim-bullseye
LABEL maintainer="glenn.s.tamkin@nasa.gov"

# --- base tools ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates bzip2 git make && \
    rm -rf /var/lib/apt/lists/*

# --- pip baseline ---
RUN python -m pip install --upgrade pip setuptools wheel
RUN pip install nvidia-pyindex nvidia-cuda-runtime-cu12

# --- micromamba (for esmpy/xesmf + compilers) ---
ARG MAMBA_VERSION=1.5.8
RUN curl -L "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}" -o /tmp/micromamba.tar.bz2 && \
    tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin bin/micromamba --strip-components=1 && \
    rm /tmp/micromamba.tar.bz2
ENV MAMBA_ROOT_PREFIX=/opt/micromamba

# Create env with things pip struggles with + compilers
RUN /usr/local/bin/micromamba create -y -n esmf -c conda-forge \
      python=3.10 esmpy xesmf libnetcdf openmpi mpi4py \
      fortran-compiler c-compiler && \
    /usr/local/bin/micromamba clean -a -y

# Make that env the default at runtime
ENV PATH=/opt/micromamba/envs/esmf/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/bin
ENV LD_LIBRARY_PATH=/opt/micromamba/envs/esmf/lib

# --- Python deps (inside the esmf env’s Python) ---
RUN python -c "import sys; print('Using Python:', sys.executable)"
RUN python -m pip install --no-cache-dir --upgrade numpy

# GraphCast / JAX / IO deps
RUN python -m pip install --no-cache-dir https://github.com/deepmind/graphcast/archive/master.zip
RUN python -m pip install --no-cache-dir --upgrade "jax[cuda12]"
RUN python -m pip install --no-cache-dir --upgrade xarray netCDF4 h5netcdf IPython

# --- Hugging Face client (for subdir download) ---
RUN python -m pip install --no-cache-dir --upgrade huggingface_hub hf_transfer
ENV HF_HUB_ENABLE_HF_TRANSFER=1

# --- Clone your repo and compile the Fortran extension ---
ARG REPO_URL="https://github.com/nasa-nccs-hpda/GenCast_FP.git"
ARG REPO_REF="deployment-operations"
WORKDIR /opt/GenCast_FP
RUN git clone --depth 1 --branch ${REPO_REF} ${REPO_URL} /opt/GenCast_FP

WORKDIR /opt/GenCast_FP/gencast_fp/preprocess
RUN python -c "import numpy as _; print('NumPy:', _.__version__)" && \
    python -m numpy.f2py -c -m eta2xprs_ eta2xprs.F90 only: eta2xprs

# Sanity check the extension
ENV PYTHONPATH=/opt/GenCast_FP
RUN python - <<'PY'
import importlib
m = importlib.import_module("gencast_fp.preprocess.eta2xprs_")
print("Built extension OK:", m.__name__)
PY

# --- Pull only the `gencast/` folder from qefm-core ---
ARG HF_QEFM_REPO="nasa-cisto-data-science-group/qefm-core"
ARG HF_QEFM_REV="main"          # optionally pin a tag or commit SHA
ARG HF_TOKEN=""                 # set if the repo is private
ENV HUGGINGFACE_HUB_TOKEN=${HF_TOKEN}

# Target location expected by the code
ENV QEFM_CORE_DIR=/opt/qefm-core
RUN mkdir -p "${QEFM_CORE_DIR}"

# Download just the subdirectory
RUN python - <<'PY'
import os
from huggingface_hub import snapshot_download

repo_id  = os.environ.get("HF_QEFM_REPO", "nasa-cisto-data-science-group/qefm-core")
revision = os.environ.get("HF_QEFM_REV", "main")
token    = os.environ.get("HUGGINGFACE_HUB_TOKEN") or None
out_dir  = os.environ["QEFM_CORE_DIR"]

print(f"[hf] Downloading {repo_id}@{revision} (only gencast/**) -> {out_dir}")
snapshot_download(
    repo_id=repo_id,
    revision=revision,
    local_dir=out_dir,
    local_dir_use_symlinks=False,
    token=token,
    allow_patterns=["gencast/**"],  # only the subdir we need
    ignore_patterns=["*.md", "*.ipynb", ".gitattributes", ".gitignore"]
)
PY

# Make the path discoverable to your app
ENV GENECAST_FP_QEFM_GENCast_DIR=${QEFM_CORE_DIR}/gencast

HEALTHCHECK NONE
ENTRYPOINT [""]
