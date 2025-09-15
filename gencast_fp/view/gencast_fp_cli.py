import sys
import time
import logging
import argparse

from gencast_fp.preprocess.fp2e5 import run_preprocess
from gencast_fp.prediction.predict_gencast import run_predict


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Process command-line args.
    parser = argparse.ArgumentParser(description="GenCast-FP Processing")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # Train
    preprocess_args = sub.add_parser("preprocess")
    preprocess_args.add_argument(
        "--start_date",
        type=str,
        required=True,
        help="Start date to process (YYYY-MM-DD)",
    )
    preprocess_args.add_argument(
        "--end_date",
        type=str,
        required=True,
        help="End date to process (YYYY-MM-DD)",
    )
    preprocess_args.add_argument(
        "--outdir",
        type=str,
        default="./output/",
        help="Output directory for the converted files",
    )
    preprocess_args.add_argument(
        "--expid",
        type=str,
        default="f5295",
        help="Experiment ID for the output files",
    )

    # Predict
    predict_args = sub.add_parser("predict")
    predict_args.add_argument(
        "--start_date",
        type=str,
        required=True,
        help="Start date to process (YYYY-MM-DD)",
    )
    predict_args.add_argument(
        "--end_date",
        type=str,
        required=True,
        help="End date to process (YYYY-MM-DD)",
    )
    predict_args.add_argument(
        "--input_dir", "-i", required=True, type=str,
        help="Preprocessed input directory")
    predict_args.add_argument(
        "--out_dir", "-o", required=True, type=str,
        help="Where to write predictions")
    predict_args.add_argument(
        "--ckpt", type=str, default=None,
        help="Path to GenCast .npz checkpoint")
    predict_args.add_argument(
        "--nsteps", type=int, default=30)
    predict_args.add_argument(
        "--res", type=float, default=1.0)
    predict_args.add_argument(
        "--ensemble", type=int, default=8)

    # Postprocess
    predict_args = sub.add_parser("postprocess")
    predict_args.add_argument("--ckpt", required=False)

    args = parser.parse_args()

    # Setup timer to monitor script execution time
    timer = time.time()

    # Execute pipeline scripts
    if args.cmd == "preprocess":
        logging.info('Starting preprocessing')
        run_preprocess(
            args.start_date,
            args.end_date,
            args.outdir,
            args.expid
        )

    elif args.cmd == "predict":
        logging.info("Starting prediction")
        out_path = run_predict(
            date=args.start_date, # TODO: this will change once we include multiday function
            input_dir=args.input_dir,
            out_dir=args.out_dir,
            ckpt_path=args.ckpt,
            res_value=args.res,
            nsteps=args.nsteps,
            ensemble_members=args.ensemble,
        )
        logging.info(f"Prediction saved: {out_path}")

    elif args.cmd == "postprocess":
        logging.info("Postprocess not implemented yet")

    else:
        raise SystemExit("Unknown command")

    logging.info(f'Took {(time.time()-timer)/60.0:.2f} min.')

    return


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
