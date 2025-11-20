# gencast_fp/predict/predict.py
import os
import logging
import dataclasses
import numpy as np
import pandas as pd
import xarray
import haiku as hk
import jax

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


def xarray_load_ds(f):
    return xarray.load_dataset(f).compute()


def load_dataset(input_dir, date, res_value, nsteps):

    # Normalize whatever we get (str, np.datetime64, Timestamp) to a Timestamp
    if not isinstance(date, pd.Timestamp):
        date = pd.to_datetime(date)

    # Match the filename convention: YYYY-MM-DDTHH
    date_str = date.strftime("%Y-%m-%dT%H")

    # Resolve dataset paths
    input_file_value = (
        f"gencast-dataset-source-geos_date-{date_str}"
        f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
    )
    input_file = os.path.join(input_dir, input_file_value)
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    return input_file


def prepare_out_dir(out_dir, date, res_value, nsteps):

    # Normalize whatever we get (str, np.datetime64, Timestamp) to a Timestamp
    if not isinstance(date, pd.Timestamp):
        date = pd.to_datetime(date)

    os.makedirs(out_dir, exist_ok=True)

    # Match the filename convention: YYYY-MM-DDTHH
    date_str = date.strftime("%Y-%m-%dT%H")

    out_file_value = (
        f"gencast-dataset-prediction-geos_date-{date_str}"
        f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
    )
    out_file = os.path.join(out_dir, out_file_value)
    return out_file


def load_ckpt_files(container_meta, ckpt_path: str = None):
    if ckpt_path is None:
        # Expect this file structure:
        #   gencast_fp/predict/predict.py (this file)
        #   ../../checkpoints/gencast/gencast-params-GenCast_1p0deg_Mini_<2019.npz
        ckpt_path = os.path.join(
            container_meta, "gencast-params-GenCast_1p0deg_Mini_<2019.npz"
        )
    diffs_file = os.path.join(
        container_meta, "gencast-stats-diffs_stddev_by_level.nc"
    )
    mean_file = os.path.join(container_meta, "gencast-stats-mean_by_level.nc")
    stddev_file = os.path.join(
        container_meta, "gencast-stats-stddev_by_level.nc"
    )
    min_file = os.path.join(container_meta, "gencast-stats-min_by_level.nc")

    with open(ckpt_path, "rb") as f:
        ckpt = checkpoint.load(f, gencast.CheckPoint)

    ckpt_and_stats = {
        "ckpt": ckpt,
        "diffs": xarray_load_ds(diffs_file),
        "mean": xarray_load_ds(mean_file),
        "stddev": xarray_load_ds(stddev_file),
        "min": xarray_load_ds(min_file),
    }

    return ckpt_and_stats


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
    res_value: float = 1.0,
    nsteps: int = 30,  # 15-day rollout (12h steps)
    ensemble_members: int = 8,
    container_meta: str = "/opt/qefm-core/gencast",
    ckpt_and_stats: dict = None,
) -> str:
    """
    Run GenCast prediction for a single date.
    Returns the path to the written NetCDF file.
    """
    logging.info("Starting prediction")
    logging.info(
        f"date={date}, input_dir={input_dir}, out_dir={out_dir}, "
        f"res={res_value}, nsteps={nsteps}, ensemble={ensemble_members}"
    )

    input_file = load_dataset(input_dir, date, res_value, nsteps)
    out_file = prepare_out_dir(out_dir, date, res_value, nsteps)

    # Skip if file already exists
    if os.path.exists(out_file):
        logging.info(f'Skipping {out_file}, prediction already exists.')
        return out_file

    # Extract model info and task info from checkpoint
    ckpt = ckpt_and_stats["ckpt"]
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

    # Load example batch to create train/eval data properly
    with open(input_file, "rb") as f:
        example_batch = xarray.load_dataset(f).compute()
    assert example_batch.dims["time"] >= 3

    # Extract train/eval tensors
    train_inputs, train_targets, train_forcings = (
        data_utils.extract_inputs_targets_forcings(
            example_batch,
            target_lead_times=slice("12h", "12h"),  # 1-AR training slice
            **dataclasses.asdict(task_config),
        )
    )
    eval_inputs, eval_targets, eval_forcings = (
        data_utils.extract_inputs_targets_forcings(
            example_batch,
            target_lead_times=slice(
                "12h", f"{(example_batch.dims['time']-2)*12}h"
            ),
            **dataclasses.asdict(task_config),
        )
    )

    logging.info(f"Devices available: {len(jax.local_devices())}")

    """
        ckpt_and_stats = {
            "ckpt": ckpt,
            "diffs": xarray_load_ds(diffs_file),
            "mean": xarray_load_ds(mean_file),
            "stddev": xarray_load_ds(stddev_file),
            "min": xarray_load_ds(min_file),
        }
    """

    # Define JIT forward pass, loss function
    def _forward_wrapped():
        predictor = _construct_wrapped_gencast(
            sampler_config,
            task_config,
            denoiser_architecture_config,
            noise_config,
            noise_encoder_config,
            ckpt_and_stats["diffs"],
            ckpt_and_stats["mean"],
            ckpt_and_stats["stddev"],
            ckpt_and_stats["min"],
        )
        return predictor

    @hk.transform_with_state
    def run_forward(inputs, targets_template, forcings):
        predictor = _forward_wrapped()
        return predictor(
            inputs, targets_template=targets_template, forcings=forcings
        )

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

    rng = jax.random.PRNGKey(0)
    rngs = np.stack(
        [jax.random.fold_in(rng, i) for i in range(ensemble_members)], axis=0
    )

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


def run_predict_multiday(
    start_date: str,
    end_date: str,
    input_dir: str,
    out_dir: str,
    res_value: float = 1.0,
    nsteps: int = 30,  # 15-day rollout (12h steps)
    ensemble_members: int = 8,
    container_meta: str = "/opt/qefm-core/gencast",
    ckpt_and_stats: dict = None,
):
    """Predict multiple days' worth of rollouts.
    Calls run_predict for each day in the start_date and end_date range."""
    #
    # start_date = np.datetime64(start_date)
    # end_date = np.datetime64(end_date)
    # date_range = np.arange(
    #     start_date, end_date + np.timedelta64(1, "D"), dtype="datetime64[D]"
    # )
    fmt = "%Y-%m-%d:%H"

    # Parse exact hour from input
    start_ts = pd.to_datetime(start_date, format=fmt)
    end_ts = pd.to_datetime(end_date,   format=fmt)

    # Generate a date range in 12-hour increments
    date_range = pd.date_range(start=start_ts, end=end_ts, freq="12h")

    for current_date in date_range:

        logging.info("======================================================")
        logging.info(f"Running prediction on date: {current_date}")
        out_fn = run_predict(
            current_date,
            input_dir,
            out_dir,
            res_value,
            nsteps,
            ensemble_members,
            container_meta,
            ckpt_and_stats,
        )
        logging.info(f"Prediction saved to file: {out_fn}")
        logging.info("======================================================")
    return out_dir


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

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    run_predict(
        date=args.date,
        input_dir=args.input_dir,
        out_dir=args.out_dir,
        ckpt_path=args.ckpt,
        res_value=args.res,
        nsteps=args.nsteps,
        ensemble_members=args.ensemble,
    )

