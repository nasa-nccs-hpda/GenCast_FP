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

    import jax
    import jax.numpy as jnp

    result = jax.jit(lambda x: x + 1)(jnp.array([1., 2.]))
    np.testing.assert_allclose(np.asarray(result), [2., 3.])
    print("Application imports and JAX CPU execution checks passed")


if __name__ == "__main__":
    main()
