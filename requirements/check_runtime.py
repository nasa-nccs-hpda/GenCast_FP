"""Exercise container runtime dependencies without model weights or a GPU."""

import importlib
from pathlib import Path
import tempfile

import numpy as np
import xarray as xr


def check_netcdf():
    """Cover explicit backends, auto-detection, file objects, and Dask reads."""
    expected = xr.Dataset(
        {"temperature": (("time", "lat", "lon"), np.arange(12.).reshape(2, 2, 3))},
        coords={
            "time": np.array(["2025-01-01", "2025-01-02"], dtype="datetime64[ns]"),
            "lat": [0., 1.],
            "lon": [10., 11., 12.],
        },
    )
    with tempfile.TemporaryDirectory() as directory:
        for engine in ("netcdf4", "h5netcdf", "scipy", None):
            path = Path(directory) / f"roundtrip-{engine}.nc"
            expected.to_netcdf(path, engine=engine)
            xr.testing.assert_allclose(expected, xr.load_dataset(path, engine=engine))
            # Match predict_gencast.py: pass an open file to automatic detection.
            with path.open("rb") as stream:
                actual = xr.load_dataset(stream).compute()
            xr.testing.assert_allclose(expected, actual)
            with xr.open_dataset(path, chunks={"time": 1}) as chunked:
                xr.testing.assert_allclose(expected, chunked.compute())
    print("NetCDF backends, file-object reads, and Dask checks passed")


def check_regridding():
    import xesmf as xe

    source = xr.Dataset(coords={"lat": [0., 1., 2.], "lon": [10., 11., 12.]})
    target = xr.Dataset(coords={"lat": [0.5, 1.5], "lon": [10.5, 11.5]})
    regridder = xe.Regridder(source, target, "bilinear")
    result = regridder(xr.DataArray(np.ones((3, 3)), dims=("lat", "lon")))
    np.testing.assert_allclose(result, 1.)
    print("ESMF/xESMF regridding check passed")


def check_ensemble_rollout():
    """Run the real sharding/pmap/rollout path with a tiny deterministic model."""
    import jax
    from graphcast import rollout, xarray_jax

    devices = jax.local_devices()
    inputs = xr.Dataset(
        {"temperature": (("time", "point"), np.array([[0., 1.], [2., 3.]], dtype=np.float32))},
        coords={"time": np.array([-12, 0], dtype="timedelta64[h]"), "point": [0, 1]},
    )
    targets = xr.Dataset(
        {"temperature": (("time", "point"), np.full((2, 2), np.nan, dtype=np.float32))},
        coords={"time": np.array([12, 24], dtype="timedelta64[h]"), "point": [0, 1]},
    )
    forcings = xr.Dataset(coords={"time": targets.time})

    def predict(rng, inputs, targets_template, forcings):
        del rng, forcings
        return targets_template.fillna(0) * 0 + inputs.isel(time=-1, drop=True) + 1

    predictor = xarray_jax.pmap(jax.jit(predict), dim="sample", devices=devices)
    samples = 2 * len(devices)
    chunks = list(rollout.chunked_prediction_generator_multiple_runs(
        predictor_fn=predictor,
        rngs=np.asarray(jax.random.split(jax.random.PRNGKey(0), samples)),
        inputs=inputs, targets_template=targets, forcings=forcings,
        num_samples=samples, pmap_devices=devices, num_steps_per_chunk=1,
    ))
    result = xr.combine_by_coords(chunks).temperature.transpose("sample", "time", "point")
    np.testing.assert_allclose(result, np.broadcast_to([[3., 4.], [4., 5.]], (samples, 2, 2)))
    np.testing.assert_array_equal(result.sample, np.arange(samples))
    np.testing.assert_array_equal(result.time, targets.time)
    print(f"Ensemble rollout passed: {samples} samples on {len(devices)} CPU devices")


def check_postprocessing():
    """Write input files, run postprocessing, and verify the CF output files."""
    import pandas as pd
    from gencast_fp.postprocess.gencast_cf import run_postprocess_multiday

    reference = pd.Timestamp("2025-01-02T12:00:00")
    for representation in ("numeric", "hours", "minutes", "timedelta", "datetime"):
        for ensemble_mean in (True, False):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                geos, predictions, output = root / "geos", root / "predictions", root / "output"
                geos.mkdir()
                predictions.mkdir()

                def fixture(hours, ensemble):
                    values = np.asarray(hours)
                    attrs = {}
                    if representation == "hours":
                        attrs = {"units": "hours", "dtype": "int64"}
                    elif representation == "minutes":
                        values = values * 60
                        attrs = {"units": "minutes"}
                    elif representation == "timedelta":
                        values = values.astype("timedelta64[h]")
                    elif representation == "datetime":
                        values = reference.to_datetime64() + values.astype("timedelta64[h]")
                    ds = xr.Dataset(
                        {"geopotential": (("batch", "time", "level", "lat", "lon"),
                                          np.full((1, 2, 2, 2, 2), 100. * 9.80665))},
                        coords={"batch": [0], "time": ("time", values, attrs),
                                "level": [500., 1000.], "lat": [1., 0.], "lon": [0., 180.]},
                    )
                    if ensemble:
                        ds = xr.concat([ds, ds * 3], dim=xr.IndexVariable("sample", [0, 1]))
                    return ds

                fixture([-12, 0], False).to_netcdf(
                    geos / "gencast-dataset-source-geos_date-2025-01-02T12_test.nc", engine="netcdf4")
                fixture([12, 24], True).to_netcdf(
                    predictions / "gencast-dataset-prediction-geos_date-2025-01-02T12_test.nc", engine="netcdf4")
                run_postprocess_multiday("2025-01-02:12", "2025-01-02:12",
                                         str(geos), str(predictions), str(output), ensemble_mean)
                files = sorted(output.rglob("*.nc"))
                assert len(files) == 4, files
                for path in files:
                    is_initial = "initial" in path.name
                    stamp = path.name.split("geos_date-")[1].split("_res-")[0]
                    with xr.open_dataset(path, engine="netcdf4") as actual:
                        np.testing.assert_array_equal(actual.time.values,
                                                      [np.datetime64(stamp, "ns")])
                        np.testing.assert_array_equal(actual.lat, [0., 1.])
                        np.testing.assert_array_equal(actual.lon, [-180., 0.])
                        np.testing.assert_array_equal(actual.lev, [1000., 500.])
                        assert "batch" not in actual.dims
                        if is_initial or ensemble_mean:
                            assert "ens" not in actual.dims
                            np.testing.assert_allclose(actual.H, 100. if is_initial else 200.)
                        else:
                            np.testing.assert_allclose(actual.H.sel(ens=0), 100.)
                            np.testing.assert_allclose(actual.H.sel(ens=1), 300.)
                    with xr.open_dataset(path, decode_cf=False) as raw:
                        np.testing.assert_array_equal(raw.time, [0.])
                        assert raw.time.dtype == np.float32
                stamps = {p.name.split("geos_date-")[1].split("_res-")[0] for p in files}
                assert stamps == {"2025-01-02T00", "2025-01-02T12", "2025-01-03T00", "2025-01-03T12"}
    print("Postprocessing passed: numeric, CF-hour/minute, timedelta, and datetime inputs; both ensemble modes")


def main():
    for module in (
        "h5py", "h5netcdf", "netCDF4", "scipy", "dask.array", "pandas",
        "haiku", "huggingface_hub", "hf_transfer", "IPython", "mpi4py.MPI",
        "gencast_fp.preprocess.fp2e5",
        "gencast_fp.prediction.predict_gencast",
        "gencast_fp.postprocess.gencast_cf",
    ):
        importlib.import_module(module)
    check_netcdf()
    check_regridding()
    check_ensemble_rollout()
    check_postprocessing()

    import jax
    import jax.numpy as jnp

    result = jax.jit(lambda x: x + 1)(jnp.array([1., 2.]))
    np.testing.assert_allclose(np.asarray(result), [2., 3.])
    print("Application imports and JAX CPU execution checks passed")


if __name__ == "__main__":
    main()
