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

# import datetime
# import math
# from typing import Optional
# import haiku as hk
import jax
from jax import random, jit, pmap
import numpy as np
import xarray
import argparse
from functools import partial

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


def parse_args():
    """Returns args parsed from command line."""
    parser = argparse.ArgumentParser(description="GenCast Mini Prediction")
    parser.add_argument(
        "--date", "-s", type=str, default="2024-12-01", help="Date to forecast"
    )
    parser.add_argument("--input_dir", "-i", type=str, help="Input directory")
    parser.add_argument("--out_dir", "-o", type=str, help="Output directory")
    return parser.parse_args()


def open_dataset(args):
    # Set up input directory and file
    date_str = args.date
    print("date_str:\n", date_str, "\n")

    # Set up input directory and file
    dataset_dir = args.input_dir
    res_value = 0.25
    steps = 30
    dataset_file_value = (
        f"gencast-dataset-prediction-geos_date-{date_str}"
        f"_res-{res_value}_levels-13_steps-{steps}.nc"
    )
    dataset_file_value = (
        f"gencast-dataset-source-geos_date-{date_str}"
        "_res-1.0_levels-13_steps-20.nc"
    )
    dataset_file = os.path.join(dataset_dir, dataset_file_value)
    if not os.path.exists(dataset_file):
        raise FileNotFoundError(f"Input file not found: {dataset_file}")

    print("dataset_file_value:\n", dataset_file_value, "\n")

    return dataset_file, dataset_file_value


def get_example_batch(dataset_file):
    def parse_file_parts(file_name):
        return dict(part.split("-", 1) for part in file_name.split("_"))

    with open(dataset_file, "rb") as f:
        example_batch = xarray.load_dataset(f).compute()

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

    return example_batch


def get_out_file(args):
    # Set up output directory
    date_str = args.date
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    res_value = 0.25
    steps = 30
    out_file_value = (
        f"gencast-dataset-prediction-geos_date-{date_str}"
        f"_res-{res_value}_levels-13_steps-{steps}.nc"
    )
    out_file = os.path.join(out_dir, out_file_value)
    return out_file


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


def load_model_ckpt(source):
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
        # SANDY -- LOAD 0.25 DEGREE CHECKPOINT
        relative_params_file = (
            "../../checkpoints/gencast/gencast-params-GenCast_0p25deg<2019.npz"
        )
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
        denoiser_architecture_config.sparse_transformer_config.mask_type = (
            "full"
        )

        print("Model description:\n", ckpt.description, "\n")
        print("Model license:\n", ckpt.license, "\n")

    return {
        "params": params,
        "state": state,
        "task_config": task_config,
        "sampler_config": sampler_config,
        "noise_config": noise_config,
        "noise_encoder_config": noise_encoder_config,
        "denoiser_architecture_config": denoiser_architecture_config,
    }


def get_train_itf(example_batch, task_config):
    return data_utils.extract_inputs_targets_forcings(
        example_batch,
        target_lead_times=slice("12h", "12h"),  # Only 1AR training.
        **dataclasses.asdict(task_config),
    )


def get_eval_itf(example_batch, task_config):
    return data_utils.extract_inputs_targets_forcings(
        example_batch,
        target_lead_times=slice(
            "12h", f"{(example_batch.dims['time']-2)*12}h"
        ),  # All but 2 input frames.
        **dataclasses.asdict(task_config),
    )


def get_norm_data():
    relative_diffs_file = (
        "../../checkpoints/gencast/gencast-stats-diffs_stddev_by_level.nc"
    )
    diffs_file = os.path.join(script_dir, relative_diffs_file)

    relative_mean_file = (
        "../../checkpoints/gencast/gencast-stats-mean_by_level.nc"
    )
    mean_file = os.path.join(script_dir, relative_mean_file)

    relative_stddev_file = (
        "../../checkpoints/gencast/gencast-stats-stddev_by_level.nc"
    )
    stddev_file = os.path.join(script_dir, relative_stddev_file)

    relative_min_file = (
        "../../checkpoints/gencast/gencast-stats-min_by_level.nc"
    )
    min_file = os.path.join(script_dir, relative_min_file)

    with open(diffs_file, "rb") as f:
        diffs_stddev_by_level = xarray.load_dataset(f).compute()
    with open(mean_file, "rb") as f:
        mean_by_level = xarray.load_dataset(f).compute()
    with open(stddev_file, "rb") as f:
        stddev_by_level = xarray.load_dataset(f).compute()
    with open(min_file, "rb") as f:
        min_by_level = xarray.load_dataset(f).compute()

    return {
        "diffs_stddev_by_level": diffs_stddev_by_level,
        "mean_by_level": mean_by_level,
        "stddev_by_level": stddev_by_level,
        "min_by_level": min_by_level,
    }


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


@jit
def run_forward(inputs, targets_template, forcings):
    predictor = construct_wrapped_gencast()
    return predictor(
        inputs, targets_template=targets_template, forcings=forcings
    )


@partial(pmap, axis_name="batch")
def run_forward_pmap(inputs, targets_template, forcings):
    return run_forward(inputs, targets_template, forcings)


@jit
def loss_fn(inputs, targets, forcings):
    predictor = construct_wrapped_gencast()
    loss, diagnostics = predictor.loss(inputs, targets, forcings)
    return xarray_tree.map_structure(
        lambda x: xarray_jax.unwrap_data(x.mean(), require_jax=True),
        (loss, diagnostics),
    )


@partial(pmap, axis_name="batch")
def loss_fn_pmap(inputs, targets, forcings):
    return loss_fn(inputs, targets, forcings)


@jit
def grads_fn(params, state, inputs, targets, forcings):
    def _aux(params, state, i, t, f):
        (loss, diagnostics), next_state = loss_fn.apply(
            params, state, random.PRNGKey(0), i, t, f
        )
        return loss, (diagnostics, next_state)

    (loss, (diagnostics, next_state)), grads = jax.value_and_grad(
        _aux, has_aux=True
    )(params, state, inputs, targets, forcings)
    return loss, diagnostics, next_state, grads


@partial(pmap, axis_name="batch")
def grads_fn_pmap(params, state, inputs, targets, forcings):
    return grads_fn(params, state, inputs, targets, forcings)


# SETUP: parse args, get dataset file, example batch, out file
args = parse_args()
dataset_file, dataset_file_value = open_dataset(args)
example_batch = get_example_batch(dataset_file)
out_file = get_out_file(args)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__name__))
print("script_dir:\n", script_dir, "\n")

# Get some random options
latent_value_options = [int(2**i) for i in range(4, 10)]
random_latent_size = 512
random_attention_type = "splash_mha"
random_mesh_size = 4
random_num_heads = 4
random_attention_k_hop = 16
random_resolution = "1p0"


# Load model ckpt and get data
source = "Checkpoint"
model_data = load_model_ckpt(source)
params = model_data["params"]
state = model_data["state"]
task_config = model_data["task_config"]
sampler_config = model_data["sampler_config"]
noise_config = model_data["noise_config"]
noise_encoder_config = model_data["noise_encoder_config"]
denoiser_architecture_config = model_data["denoiser_architecture_config"]

# Get inputs, targets, forcings for train/eval splits
train_inputs, train_targets, train_forcings = get_train_itf(
    example_batch, task_config
)
eval_inputs, eval_targets, eval_forcings = get_eval_itf(
    example_batch, task_config
)

print("All Examples:  ", example_batch.dims.mapping)
print("Train Inputs:  ", train_inputs.dims.mapping)
print("Train Targets: ", train_targets.dims.mapping)
print("Train Forcings:", train_forcings.dims.mapping)
print("Eval Inputs:   ", eval_inputs.dims.mapping)
print("Eval Targets:  ", eval_targets.dims.mapping)
print("Eval Forcings: ", eval_forcings.dims.mapping)

# Get normalization values
norm_data = get_norm_data()
diffs_stddev_by_level = norm_data["diffs_stddev_by_level"]
mean_by_level = norm_data["mean_by_level"]
stddev_by_level = norm_data["stddev_by_level"]
min_by_level = norm_data["min_by_level"]

# Initialize distributed processing
if "SLURM_PROCID" in os.environ:
    coordinator = f"{os.environ['MASTER_ADDR']}:{os.environ['MASTER_PORT']}"
    jax.distributed.initialize(
        coordinator_address=coordinator,
        num_processes=int(os.environ["WORLD_SIZE"]),
        process_id=int(os.environ["SLURM_PROCID"]),
    )

print(
    f"Process {jax.process_index()} ready on "
    f"{jax.local_device_count()} local devices."
)

rng = random.PRNGKey(0)
num_ensemble_members = 8
rngs = np.stack(
    [random.fold_in(rng, i) for i in range(num_ensemble_members)], axis=0
)

eval_inputs = jax.device_put_replicated(eval_inputs, jax.local_devices())
eval_targets = jax.device_put_replicated(eval_targets, jax.local_devices())
eval_forcings = jax.device_put_replicated(eval_forcings, jax.local_devices())

# Run chunked prediction using pmapped functions
chunks = []

chunked_pred_generator = rollout.chunked_prediction_generator_multiple_runs(
    # Use pmapped version to parallelise across devices.
    predictor_fn=run_forward_pmap,
    rngs=rngs,
    inputs=eval_inputs,
    targets_template=eval_targets * np.nan,
    forcings=eval_forcings,
    num_steps_per_chunk=1,
    num_samples=num_ensemble_members,
    pmap_devices=jax.local_devices(),
)

for chunk in chunked_pred_generator:
    chunks.append(chunk)
predictions = xarray.combine_by_coords(chunks)
predictions.to_netcdf(out_file)
print("Predictions for 15 days complete. Out file:\n", out_file, "\n")
