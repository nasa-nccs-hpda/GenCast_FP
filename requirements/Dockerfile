FROM python:3.10-slim-bullseye
LABEL maintainer="glenn.s.tamkin@nasa.gov"

# system deps you were already using (keep as-is)
RUN pip install --upgrade setuptools pip wheel
RUN pip install nvidia-pyindex
RUN pip install nvidia-cuda-runtime-cu12

# ---- add micromamba (minimal) ----
ARG MAMBA_VERSION=1.5.8
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates bzip2 && \
    rm -rf /var/lib/apt/lists/*
RUN curl -L "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}" -o /tmp/micromamba.tar.bz2 && \
    tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin bin/micromamba --strip-components=1 && \
    rm /tmp/micromamba.tar.bz2
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
# Create a lightweight env with only what pip can’t provide
RUN /usr/local/bin/micromamba create -y -n esmf -c conda-forge \
      python=3.10 esmpy xesmf libnetcdf openmpi mpi4py && \
    /usr/local/bin/micromamba clean -a -y
# Use the conda env’s python/pkgs at runtime
ENV PATH=$MAMBA_ROOT_PREFIX/envs/esmf/bin:$PATH
ENV LD_LIBRARY_PATH=$MAMBA_ROOT_PREFIX/envs/esmf/lib:$LD_LIBRARY_PATH

# ---- your existing pip installs (these now run inside the esmf env’s python) ----
# GraphCast dependencies
RUN pip install --no-cache-dir https://github.com/deepmind/graphcast/archive/master.zip

# GPU deps (JAX CUDA)
RUN pip install --no-cache-dir --upgrade "jax[cuda12]"

# GenCast deps
RUN pip install --no-cache-dir --upgrade xarray netCDF4 h5netcdf IPython

HEALTHCHECK NONE
ENTRYPOINT [""]
