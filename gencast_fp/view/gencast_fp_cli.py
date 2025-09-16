import os
import sys
import time
import logging
import argparse
from pathlib import Path

from gencast_fp.preprocess.fp2e5 import run_preprocess
from gencast_fp.prediction.predict_gencast import (
    run_predict_multiday,
    load_ckpt_files,
)
from gencast_fp.postprocess.gencast_cf import run_postprocess_multiday

# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(description="GenCast-FP Processing")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # ---------- preprocess ----------
    preprocess_args = sub.add_parser("preprocess", help="Run preprocessing only")
    preprocess_args.add_argument("--start_date", type=str, required=True, help="YYYY-MM-DD")
    preprocess_args.add_argument("--end_date",   type=str, required=True, help="YYYY-MM-DD")
    preprocess_args.add_argument("--output_dir",     type=str, default="./output/preprocess",
                                 help="Output directory for preprocessed files")
    preprocess_args.add_argument("--expid",      type=str, default="f5295",
                                 help="Experiment ID for the output files")

    # ---------- predict ----------
    predict_args = sub.add_parser("predict", help="Run prediction only")
    predict_args.add_argument("--start_date", type=str, required=True, help="YYYY-MM-DD")
    predict_args.add_argument("--end_date",   type=str, required=True, help="YYYY-MM-DD")
    predict_args.add_argument("--input_dir",  "-i", required=True, type=str,
                              help="Preprocessed input directory")
    predict_args.add_argument("--output_dir",    "-o", required=True, type=str,
                              help="Where to write predictions")
    predict_args.add_argument("--ckpt",       type=str, default=None,
                              help="Path to GenCast .npz checkpoint (overrides container default)")
    predict_args.add_argument("--nsteps",     type=int, default=30)
    predict_args.add_argument("--res",        type=float, default=1.0)
    predict_args.add_argument("--ensemble",   type=int, default=8)
    predict_args.add_argument("--container_meta",  type=str, default="/opt/qefm-core/gencast",
                          help="Where to load default ckpt/configs if --ckpt not passed")

    # ---------- postprocess ----------
    post_args = sub.add_parser("postprocess", help="Run postprocessing only")
    post_args.add_argument("--start_date", type=str, required=True, help="Start date (YYYY-MM-DD)")
    post_args.add_argument("--end_date",   type=str, required=True, help="End date (YYYY-MM-DD)")
    post_args.add_argument("--input_dir",   type=str, required=True,
                        help="Directory with GEOS inputs (for initial conditions)")
    post_args.add_argument("--predictions_dir",   type=str, required=True,
                        help="Directory with GenCast predictions")
    post_args.add_argument("--output_dir", type=str, default="./output/postprocess",
                        help="Directory for CF-compliant NetCDF outputs")
    post_args.add_argument("--no_ens_mean", action="store_true",
                        help="Disable ensemble mean (keep all ensemble members)")

    # ---------- run (all-in-one) ----------
    run_args = sub.add_parser("run", help="Run preprocess → predict → postprocess")
    run_args.add_argument("--start_date", type=str, required=True, help="YYYY-MM-DD")
    run_args.add_argument("--end_date",   type=str, required=True, help="YYYY-MM-DD")

    run_args.add_argument("--output_dir", type=str, default="./output/preprocess",
                          help="Directory for preprocess outputs (becomes predict input)")
    run_args.add_argument("--expid",          type=str, default="f5295",
                          help="Experiment ID used during preprocessing")

    run_args.add_argument("--ckpt",     type=str, default=None,
                          help="Path to GenCast .npz checkpoint (overrides container default)")
    run_args.add_argument("--nsteps",   type=int, default=30)
    run_args.add_argument("--res",      type=float, default=1.0)
    run_args.add_argument("--ensemble", type=int, default=8)

    run_args.add_argument("--skip_preprocess", action="store_true",
                          help="Skip preprocess (assumes --preprocess_dir already exists)")
    run_args.add_argument("--skip_predict",    action="store_true",
                          help="Skip predict")
    run_args.add_argument("--skip_post",       action="store_true",
                          help="Skip postprocess")
    run_args.add_argument("--container_meta",  type=str, default="/opt/qefm-core/gencast",
                          help="Where to load default ckpt/configs if --ckpt not passed")

    args = parser.parse_args()
    t0 = time.time()

    if args.cmd == "preprocess":
        logging.info("Starting preprocessing")

        run_preprocess(
            args.start_date, args.end_date, args.out_dir, args.expid)

    elif args.cmd == "predict":

        logging.info("Starting prediction...")

        if args.ckpt:
            ckpt_and_stats = {"ckpt": args.ckpt}
        else:
            ckpt_and_stats = load_ckpt_files("/opt/qefm-core/gencast")

        out_path = run_predict_multiday(
            start_date=args.start_date,
            end_date=args.end_date,
            input_dir=args.input_dir,
            out_dir=args.out_dir,
            res_value=args.res,
            nsteps=args.nsteps,
            ensemble_members=args.ensemble,
            container_meta=args.container_meta,
            ckpt_and_stats=ckpt_and_stats,
        )
        logging.info(f"Prediction saved: {out_path}")

    elif args.cmd == "postprocess":

        run_postprocess_multiday(
            start_date=args.start_date,
            end_date=args.end_date,
            geos_dir=args.input_dir,
            pred_dir=args.predictions_dir,
            post_out_dir=args.output_dir,
            ens_mean=args.ens_mean,
        )

    elif args.cmd == "run":

        # Setting up directories
        preprocess_output_dir = os.path.join(args.output_dir, '1-preprocessed')
        prediction_output_dir = os.path.join(args.output_dir, '2-predictions')
        postprocess_output_dir = os.path.join(args.output_dir, '3-postprocessed')
        for edir in [
                    preprocess_output_dir, prediction_output_dir,
                    postprocess_output_dir
                ]:
            os.makedirs(edir, exist_ok=True)

        # 1) Preprocess
        if not args.skip_preprocess:
            logging.info(f"[1/3] Preprocess → {preprocess_output_dir}")
            run_preprocess(
                args.start_date, args.end_date,
                preprocess_output_dir, args.expid
            )
        else:
            logging.info("[1/3] Skipping preprocess")

        # 2) Predict
        if not args.skip_predict:
            logging.info(f"[2/3] Prediction → {prediction_output_dir}")
            if args.ckpt:
                ckpt_and_stats = {"ckpt": args.ckpt}
            else:
                ckpt_and_stats = load_ckpt_files(args.container_meta)

            out_path = run_predict_multiday(
                start_date=args.start_date,
                end_date=args.end_date,
                input_dir=preprocess_output_dir,
                out_dir=prediction_output_dir,
                res_value=args.res,
                nsteps=args.nsteps,
                ensemble_members=args.ensemble,
                container_meta=args.container_meta,
                ckpt_and_stats=ckpt_and_stats,
            )
            logging.info(f"Prediction saved: {out_path}")
        else:
            logging.info("[2/3] Skipping predict")

        # 3) Postprocess
        if not args.skip_post:

            logging.info(f"[3/3] Postprocess → {postprocess_output_dir}")

            run_postprocess_multiday(
                start_date=args.start_date,
                end_date=args.end_date,
                geos_dir=preprocess_output_dir,
                pred_dir=prediction_output_dir,
                post_out_dir=postprocess_output_dir,
                ens_mean=args.ens_mean,
            )

        else:
            logging.info("[3/3] Skipping postprocess")

    else:
        raise SystemExit("Unknown command")

    logging.info(f"Took {(time.time()-t0)/60.0:.2f} min.")


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
