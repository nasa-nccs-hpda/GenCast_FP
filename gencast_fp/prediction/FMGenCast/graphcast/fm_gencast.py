#!/usr/bin/env python
# coding: utf-8

# > Copyright 2024 DeepMind Technologies Limited.
# >
# > Licensed under the Apache License, Version 2.0 (the "License");
# > you may not use this file except in compliance with the License.
# > You may obtain a copy of the License at
# >
# >      http://www.apache.org/licenses/LICENSE-2.0
# >
# > Unless required by applicable law or agreed to in writing, software
# > distributed under the License is distributed on an "AS-IS" BASIS,
# > WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# > See the License for the specific language governing permissions and
# > limitations under the License.

# @title Imports

import dataclasses
import datetime
import math
from typing import Optional
import haiku as hk
import jax
import numpy as np
import xarray
import argparse

from graphcast import rollout
from graphcast import xarray_jax
from graphcast import normalization
from graphcast import checkpoint
from graphcast import data_utils
from graphcast import xarray_tree
from graphcast import gencast
from graphcast import denoiser
from graphcast import nan_cleaning
import os

parser = argparse.ArgumentParser(description="GenCast Mini Prediction")
parser.add_argument(
    "--date", "-s", type=str, default="2024-12-01", help="Date to forecast"
)
parser.add_argument("--input_dir", "-i", type=str, help="Input directory")
parser.add_argument("--out_dir", "-o", type=str, help="Output directory")


args = parser.parse_args()
# Set up input directory and file
date_str = args.date
print("date_str:\n", date_str, "\n")


# Set up input directory and file
dataset_dir = args.input_dir
res_value = 1.0  # 1.0 resolution
nsteps = 30  # 15 day rollout
dataset_file_value = (
    f"gencast-dataset-source-geos_date-{date_str}"
    f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
)
dataset_file = os.path.join(dataset_dir, dataset_file_value)
if not os.path.exists(dataset_file):
    raise FileNotFoundError(f"Input file not found: {dataset_file}")

# Set up output directory
out_dir = args.out_dir
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file_value = (
    f"gencast-dataset-prediction-geos_date-{date_str}"
    f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
)
out_file = os.path.join(out_dir, out_file_value)


# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__name__))
print("script_dir:\n", script_dir, "\n")


latent_value_options = [int(2**i) for i in range(4, 10)]

random_latent_size = 512
random_attention_type = "splash_mha"
random_mesh_size = 4
random_num_heads = 4
random_attention_k_hop = 16
random_resolution = "1p0"


def update_latent_options(*args):
    def _latent_valid_for_attn(attn, latent, heads):
        head_dim, rem = divmod(latent, heads)
        if rem != 0:
            return False
        # Required for splash attn.
        if head_dim % 128 != 0:
            return attn != "splash_mha"
        return True

    attn = random_attention_type.value
    heads = random_num_heads.value
    random_latent_size.options = [
        latent
        for latent in latent_value_options
        if _latent_valid_for_attn(attn, latent, heads)
    ]


# @title Load the model
source = "Checkpoint"
if source == "Random":
    params = None  # Filled in below
    state = {}
    task_config = gencast.TASK
    # Use default values.
    sampler_config = gencast.SamplerConfig()
    noise_config = gencast.NoiseConfig()
    noise_encoder_config = denoiser.NoiseEncoderConfig()
    # Configure, otherwise use default values.
    denoiser_architecture_config = denoiser.DenoiserArchitectureConfig(
        sparse_transformer_config=denoiser.SparseTransformerConfig(
            attention_k_hop=random_attention_k_hop.value,
            attention_type=random_attention_type.value,
            d_model=random_latent_size.value,
            num_heads=random_num_heads.value,
        ),
        mesh_size=random_mesh_size.value,
        latent_size=random_latent_size.value,
    )
else:
    assert source == "Checkpoint"
    params_file_value = "GenCast 1p0deg Mini <2019.npz"
    relative_params_file = "../../checkpoints/gencast/gencast-params-GenCast_1p0deg_Mini_<2019.npz"
    absolute_path = os.path.join(script_dir, relative_params_file)
    print("absolute_path:\n", absolute_path, "\n")
    params_file = absolute_path
    with open(params_file, "rb") as f:
        print(params_file)
        ckpt = checkpoint.load(f, gencast.CheckPoint)
    params = ckpt.params
    state = {}

    task_config = ckpt.task_config
    sampler_config = ckpt.sampler_config
    noise_config = ckpt.noise_config
    noise_encoder_config = ckpt.noise_encoder_config
    denoiser_architecture_config = ckpt.denoiser_architecture_config

    denoiser_architecture_config.sparse_transformer_config.attention_type = (
        "triblockdiag_mha"
    )
    denoiser_architecture_config.sparse_transformer_config.mask_type = "full"

    print("Model description:\n", ckpt.description, "\n")
    print("Model license:\n", ckpt.license, "\n")


# ## Load the example data
def parse_file_parts(file_name):
    return dict(part.split("-", 1) for part in file_name.split("_"))


def data_valid_for_model(file_name: str, params_file_name: str):
    """Check data type and resolution matches."""
    data_file_parts = parse_file_parts(file_name.removesuffix(".nc"))
    data_res = data_file_parts["res"].replace(".", "p")
    if source == "Random":
        return random_resolution.value == data_res
    res_matches = data_res in params_file_name.lower()
    source_matches = "Operational" in params_file_name
    if data_file_parts["source"] == "era5":
        source_matches = not source_matches
    return res_matches and source_matches


# @title Load weather data
# dataset_dir = "/discover/nobackup/jli30/GenCast_FP/output_test"
# dataset_file_value = f"gencast-dataset-source-geos_date-2024-12-01_res-1.0_levels-13_steps-20.nc"
# dataset_file = os.path.join(dataset_dir, dataset_file_value)
print("dataset_file_value:\n", dataset_file_value, "\n")
# with gcs_bucket.blob(dir_prefix + f"dataset/{dataset_file_value}").open("rb") as f:
with open(dataset_file, "rb") as f:
    example_batch = xarray.load_dataset(f).compute()
##example_batch = xarray.open_dataset(dataset_file)

assert example_batch.dims["time"] >= 3  # 2 for input, >=1 for targets

print(
    ", ".join(
        [
            f"{k}: {v}"
            for k, v in parse_file_parts(
                dataset_file_value.removesuffix(".nc")
            ).items()
        ]
    )
)

# print(example_batch['2m_temperature'].isel(time=0).squeeze().to_numpy())
##example_batch
# @title Extract training and eval data

train_inputs, train_targets, train_forcings = (
    data_utils.extract_inputs_targets_forcings(
        example_batch,
        target_lead_times=slice("12h", "12h"),  # Only 1AR training.
        **dataclasses.asdict(task_config),
    )
)

eval_inputs, eval_targets, eval_forcings = (
    data_utils.extract_inputs_targets_forcings(
        example_batch,
        target_lead_times=slice(
            "12h", f"{(example_batch.dims['time']-2)*12}h"
        ),  # All but 2 input frames.
        **dataclasses.asdict(task_config),
    )
)

print("All Examples:  ", example_batch.dims.mapping)
print("Train Inputs:  ", train_inputs.dims.mapping)
print("Train Targets: ", train_targets.dims.mapping)
print("Train Forcings:", train_forcings.dims.mapping)
print("Eval Inputs:   ", eval_inputs.dims.mapping)
print("Eval Targets:  ", eval_targets.dims.mapping)
print("Eval Forcings: ", eval_forcings.dims.mapping)

# @title Load normalization data
relative_diffs_file = (
    "../../checkpoints/gencast/gencast-stats-diffs_stddev_by_level.nc"
)
diffs_file = os.path.join(script_dir, relative_diffs_file)

relative_mean_file = "../../checkpoints/gencast/gencast-stats-mean_by_level.nc"
mean_file = os.path.join(script_dir, relative_mean_file)

relative_stddev_file = (
    "../../checkpoints/gencast/gencast-stats-stddev_by_level.nc"
)
stddev_file = os.path.join(script_dir, relative_stddev_file)

relative_min_file = "../../checkpoints/gencast/gencast-stats-min_by_level.nc"
min_file = os.path.join(script_dir, relative_min_file)

with open(diffs_file, "rb") as f:
    diffs_stddev_by_level = xarray.load_dataset(f).compute()
with open(mean_file, "rb") as f:
    mean_by_level = xarray.load_dataset(f).compute()
with open(stddev_file, "rb") as f:
    stddev_by_level = xarray.load_dataset(f).compute()
with open(min_file, "rb") as f:
    min_by_level = xarray.load_dataset(f).compute()


def construct_wrapped_gencast():
    """Constructs and wraps the GenCast Predictor."""
    predictor = gencast.GenCast(
        sampler_config=sampler_config,
        task_config=task_config,
        denoiser_architecture_config=denoiser_architecture_config,
        noise_config=noise_config,
        noise_encoder_config=noise_encoder_config,
    )

    predictor = normalization.InputsAndResiduals(
        predictor,
        diffs_stddev_by_level=diffs_stddev_by_level,
        mean_by_level=mean_by_level,
        stddev_by_level=stddev_by_level,
    )

    predictor = nan_cleaning.NaNCleaner(
        predictor=predictor,
        reintroduce_nans=True,
        fill_value=min_by_level,
        var_to_clean="sea_surface_temperature",
    )

    return predictor


@hk.transform_with_state
def run_forward(inputs, targets_template, forcings):
    predictor = construct_wrapped_gencast()
    return predictor(
        inputs, targets_template=targets_template, forcings=forcings
    )


@hk.transform_with_state
def loss_fn(inputs, targets, forcings):
    predictor = construct_wrapped_gencast()
    loss, diagnostics = predictor.loss(inputs, targets, forcings)
    return xarray_tree.map_structure(
        lambda x: xarray_jax.unwrap_data(x.mean(), require_jax=True),
        (loss, diagnostics),
    )


def grads_fn(params, state, inputs, targets, forcings):
    def _aux(params, state, i, t, f):
        (loss, diagnostics), next_state = loss_fn.apply(
            params, state, jax.random.PRNGKey(0), i, t, f
        )
        return loss, (diagnostics, next_state)

    (loss, (diagnostics, next_state)), grads = jax.value_and_grad(
        _aux, has_aux=True
    )(params, state, inputs, targets, forcings)
    return loss, diagnostics, next_state, grads


if params is None:
    init_jitted = jax.jit(loss_fn.init)
    params, state = init_jitted(
        rng=jax.random.PRNGKey(0),
        inputs=train_inputs,
        targets=train_targets,
        forcings=train_forcings,
    )
# GST orig below
grads_fn_jitted = jax.jit(grads_fn)
run_forward_jitted = jax.jit(
    lambda rng, i, t, f: run_forward.apply(params, state, rng, i, t, f)[0]
)
# We also produce a pmapped version for running in parallel.
run_forward_pmap = xarray_jax.pmap(run_forward_jitted, dim="sample")


print(f"Number of local devices {len(jax.local_devices())}")

# @title Autoregressive rollout (loop in python)

print("Inputs:  ", eval_inputs.dims.mapping)
print("Targets: ", eval_targets.dims.mapping)
print("Forcings:", eval_forcings.dims.mapping)

num_ensemble_members = 8  # @param int
rng = jax.random.PRNGKey(0)
# We fold-in the ensemble member, this way the first N members should always
# match across different runs which use take the same inputs, regardless of
# total ensemble size.
rngs = np.stack(
    [jax.random.fold_in(rng, i) for i in range(num_ensemble_members)], axis=0
)

chunks = []

# GST
for chunk in rollout.chunked_prediction_generator_multiple_runs(
    # Use pmapped version to parallelise across devices.
    predictor_fn=run_forward_pmap,
    rngs=rngs,
    inputs=eval_inputs,
    targets_template=eval_targets * np.nan,
    forcings=eval_forcings,
    num_steps_per_chunk=1,
    num_samples=num_ensemble_members,
    pmap_devices=jax.local_devices(),
):
    chunks.append(chunk)
predictions = xarray.combine_by_coords(chunks)
predictions.to_netcdf(out_file)
print("Predictions computed for 10 days out_file:\n", out_file, "\n")
