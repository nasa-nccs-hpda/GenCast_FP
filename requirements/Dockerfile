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

# --- Python deps (these run inside the esmf env’s python) ---
RUN python -c "import sys; print('Using Python:', sys.executable)"
# Ensure numpy is present for f2py
RUN python -m pip install --no-cache-dir --upgrade numpy

# GraphCast / JAX / IO deps
RUN python -m pip install --no-cache-dir https://github.com/deepmind/graphcast/archive/master.zip
RUN python -m pip install --no-cache-dir --upgrade "jax[cuda12]"
RUN python -m pip install --no-cache-dir --upgrade xarray netCDF4 h5netcdf IPython

# --- Clone your repo and compile the Fortran extension ---
# Customize these build args to your repo
ARG REPO_URL="https://github.com/nasa-nccs-hpda/GenCast_FP.git"
ARG REPO_REF="deployment-operations"  # or a specific commit/branch/tag
WORKDIR /opt/GenCast_FP
RUN git clone --depth 1 --branch ${REPO_REF} ${REPO_URL} /opt/GenCast_FP

# Compile eta2xprs.F90 into eta2xprs_.so for CPython 3.10
WORKDIR /opt/GenCast_FP/gencast_fp/preprocess
# f2py will place eta2xprs_.cpython-310-...so in this directory
RUN python -c "import numpy as _; print('NumPy:', _.__version__)" && \
    python -m numpy.f2py -c -m eta2xprs_ eta2xprs.F90 only: eta2xprs

# Sanity: import the compiled module at build time
ENV PYTHONPATH=/opt/GenCast_FP
RUN python - <<'PY'
import importlib
m = importlib.import_module("gencast_fp.preprocess.eta2xprs_")
print("Built extension OK:", m.__name__)
PY

# (optional) place checkpoints/stats and point QEFM_CORE_DIR there
# ENV QEFM_CORE_DIR=/opt/qefm-core
# ... add a snapshot_download or COPY step if you want to bake assets in

HEALTHCHECK NONE
ENTRYPOINT [""]
