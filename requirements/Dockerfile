## Start with Debian slim-bullseye - https://medium.com/@arif.rahman.rhm/choosing-the-right-python-docker-image-slim-buster-vs-alpine-vs-slim-bullseye-5586bac8b4c9
# Choose Slim Bullseye if:
# You want a newer base image with reasonable size.
# You prefer the familiarity of Debian while benefiting from some size reduction.

FROM python:3.10-slim-bullseye
LABEL maintainer="glenn.s.tamkin@nasa.gov" 

RUN pip install --upgrade setuptools pip wheel
RUN pip install nvidia-pyindex
RUN pip install nvidia-cuda-runtime-cu12

# GraphCast dependencies
RUN pip install --no-cache-dir https://github.com/deepmind/graphcast/archive/master.zip

# GPU dependencies
RUN pip install --upgrade "jax[cuda12]"

# GenCast dependencies
RUN pip install --upgrade xarray netCDF4 h5netcdf

# Regridding dependencies for preprocess
RUN apt-get update && apt-get install -y --no-install-recommends \
        libnetcdf-dev \
        libopenmpi-dev \
        gfortran \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir esmpy xesmf

RUN pip install IPython
#RUN pip install --upgrade bitsandbytes

# Aurora dependencies
#RUN pip install cdsapi

# Just cuz....https://docs.jax.dev/en/latest/installation.html#install-nvidia-gpu
#RUN pip install --upgrade "jax[cuda12]"

HEALTHCHECK NONE
ENTRYPOINT [""]