import sys
import time
import logging
import argparse

from gencast_fp.preprocess.fp2e5 import run_preprocess


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
    predict_args.add_argument("--ckpt", required=False)
    predict_args.add_argument("--wfs", required=False)
    predict_args.add_argument("--sci", default=None)
    predict_args.add_argument("--out-dir", default="pred_out")
    predict_args.add_argument("--batch-size", type=int, default=128)
    predict_args.add_argument("--num-workers", type=int, default=4)
    predict_args.add_argument("--seed", type=int, default=42)
    predict_args.add_argument("--devices", type=int, default=1)
    predict_args.add_argument("--precision", type=int, default=32)
    predict_args.add_argument("--save-pngs", type=int, default=0)
    predict_args.add_argument("--max-pngs", type=int, default=0)

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

    #elif args.cmd == "predict":
    #    pipeline.predict(args)
    #elif args.cmd == "detect":
    #    pipeline.detect(args)
    #else:
    #    raise SystemExit("Unknown command")

    logging.info(f'Took {(time.time()-timer)/60.0:.2f} min.')

    return


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
