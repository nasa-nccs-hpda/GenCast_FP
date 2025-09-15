# gencast_fp/predict/predict.py
import os
import dataclasses
import logging
import numpy as np
import xarray
import haiku as hk
import jax

from typing import Optional

from graphcast import rollout
from graphcast import xarray_jax
from graphcast import normalization
from graphcast import checkpoint
from graphcast import data_utils
from graphcast import xarray_tree
from graphcast import gencast
from graphcast import denoiser
from graphcast import nan_cleaning


def _parse_file_parts(file_name: str):
    return dict(part.split("-", 1) for part in file_name.split("_"))


def _construct_wrapped_gencast(
    sampler_config,
    task_config,
    denoiser_architecture_config,
    noise_config,
    noise_encoder_config,
    diffs_stddev_by_level,
    mean_by_level,
    stddev_by_level,
    min_by_level,
):
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


def run_predict(
    date: str,
    input_dir: str,
    out_dir: str,
    ckpt_path: Optional[str] = None,
    res_value: float = 1.0,
    nsteps: int = 30,                # 15-day rollout (12h steps)
    ensemble_members: int = 8,
) -> str:
    """
    Run GenCast prediction for a single date.

    Returns the path to the written NetCDF file.
    """
    logging.info("Starting prediction")
    logging.info(f"date={date}, input_dir={input_dir}, out_dir={out_dir}, "
                 f"res={res_value}, nsteps={nsteps}, ensemble={ensemble_members}")

    # Resolve dataset paths
    dataset_file_value = (
        f"gencast-dataset-source-geos_date-{date}"
        f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
    )
    dataset_file = os.path.join(input_dir, dataset_file_value)
    if not os.path.exists(dataset_file):
        raise FileNotFoundError(f"Input file not found: {dataset_file}")
    os.makedirs(out_dir, exist_ok=True)
    out_file_value = (
        f"gencast-dataset-prediction-geos_date-{date}"
        f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
    )
    out_file = os.path.join(out_dir, out_file_value)

    # Resolve checkpoints and stats
    # Default to repo-relative checkpoints if ckpt_path not given.
    if ckpt_path is None:
        # Expect this file structure:
        #   gencast_fp/predict/predict.py (this file)
        #   ../../checkpoints/gencast/gencast-params-GenCast_1p0deg_Mini_<2019.npz
        script_dir = os.path.dirname(os.path.abspath(__file__))
        ckpt_path = os.path.join(
            script_dir,
            "../../checkpoints/gencast/gencast-params-GenCast_1p0deg_Mini_<2019.npz",
        )
    diffs_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../checkpoints/gencast/gencast-stats-diffs_stddev_by_level.nc",
    )
    mean_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../checkpoints/gencast/gencast-stats-mean_by_level.nc",
    )
    stddev_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../checkpoints/gencast/gencast-stats-stddev_by_level.nc",
    )
    min_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../checkpoints/gencast/gencast-stats-min_by_level.nc",
    )

    logging.info(f"Checkpoint: {ckpt_path}")
    logging.info("Loading checkpoint and configs...")
    with open(ckpt_path, "rb") as f:
        ckpt = checkpoint.load(f, gencast.CheckPoint)
    params = ckpt.params
    state = {}
    task_config = ckpt.task_config
    sampler_config = ckpt.sampler_config
    noise_config = ckpt.noise_config
    noise_encoder_config = ckpt.noise_encoder_config
    denoiser_architecture_config = ckpt.denoiser_architecture_config
    # Keep the same overrides you had:
    denoiser_architecture_config.sparse_transformer_config.attention_type = "triblockdiag_mha"
    denoiser_architecture_config.sparse_transformer_config.mask_type = "full"

    logging.info("Loading dataset and normalization stats...")
    with open(dataset_file, "rb") as f:
        example_batch = xarray.load_dataset(f).compute()
    assert example_batch.dims["time"] >= 3

    with open(diffs_file, "rb") as f:
        diffs_stddev_by_level = xarray.load_dataset(f).compute()
    with open(mean_file, "rb") as f:
        mean_by_level = xarray.load_dataset(f).compute()
    with open(stddev_file, "rb") as f:
        stddev_by_level = xarray.load_dataset(f).compute()
    with open(min_file, "rb") as f:
        min_by_level = xarray.load_dataset(f).compute()

    # Extract train/eval tensors
    train_inputs, train_targets, train_forcings = data_utils.extract_inputs_targets_forcings(
        example_batch,
        target_lead_times=slice("12h", "12h"),  # 1-AR training slice
        **dataclasses.asdict(task_config),
    )
    eval_inputs, eval_targets, eval_forcings = data_utils.extract_inputs_targets_forcings(
        example_batch,
        target_lead_times=slice("12h", f"{(example_batch.dims['time']-2)*12}h"),
        **dataclasses.asdict(task_config),
    )

    logging.info(f"Devices available: {len(jax.local_devices())}")

    def _forward_wrapped():
        predictor = _construct_wrapped_gencast(
            sampler_config,
            task_config,
            denoiser_architecture_config,
            noise_config,
            noise_encoder_config,
            diffs_stddev_by_level,
            mean_by_level,
            stddev_by_level,
            min_by_level,
        )
        return predictor

    @hk.transform_with_state
    def run_forward(inputs, targets_template, forcings):
        predictor = _forward_wrapped()
        return predictor(inputs, targets_template=targets_template, forcings=forcings)

    @hk.transform_with_state
    def loss_fn(inputs, targets, forcings):
        predictor = _forward_wrapped()
        loss, diagnostics = predictor.loss(inputs, targets, forcings)
        return xarray_tree.map_structure(
            lambda x: xarray_jax.unwrap_data(x.mean(), require_jax=True),
            (loss, diagnostics),
        )

    if params is None:
        init_jitted = jax.jit(loss_fn.init)
        params, state = init_jitted(
            rng=jax.random.PRNGKey(0),
            inputs=train_inputs,
            targets=train_targets,
            forcings=train_forcings,
        )

    run_forward_jitted = jax.jit(
        lambda rng, i, t, f: run_forward.apply(params, state, rng, i, t, f)[0]
    )
    run_forward_pmap = xarray_jax.pmap(run_forward_jitted, dim="sample")

    # Rollout
    rng = jax.random.PRNGKey(0)
    rngs = np.stack([jax.random.fold_in(rng, i) for i in range(ensemble_members)], axis=0)

    logging.info("Starting autoregressive rollout...")
    chunks = []
    for chunk in rollout.chunked_prediction_generator_multiple_runs(
        predictor_fn=run_forward_pmap,
        rngs=rngs,
        inputs=eval_inputs,
        targets_template=eval_targets * np.nan,
        forcings=eval_forcings,
        num_steps_per_chunk=1,
        num_samples=ensemble_members,
        pmap_devices=jax.local_devices(),
    ):
        chunks.append(chunk)
    predictions = xarray.combine_by_coords(chunks)
    predictions.to_netcdf(out_file)

    logging.info(f"Prediction written: {out_file}")
    return out_file


if __name__ == "__main__":
    # Optional standalone CLI for this module
    import argparse
    parser = argparse.ArgumentParser(description="GenCast Mini Prediction")
    parser.add_argument("--date", "-s", type=str, required=True)
    parser.add_argument("--input_dir", "-i", type=str, required=True)
    parser.add_argument("--out_dir", "-o", type=str, required=True)
    parser.add_argument("--ckpt", type=str, default=None)
    parser.add_argument("--nsteps", type=int, default=30)
    parser.add_argument("--res", type=float, default=1.0)
    parser.add_argument("--ensemble", type=int, default=8)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    run_predict(
        date=args.date,
        input_dir=args.input_dir,
        out_dir=args.out_dir,
        ckpt_path=args.ckpt,
        res_value=args.res,
        nsteps=args.nsteps,
        ensemble_members=args.ensemble,
    )
