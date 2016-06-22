#! /usr/bin/env python
import logging
import os

import Logging
import PB8KFileInfo

BETA="BETA"
ZSCORE="ZSCORE"

BEST="BEST"
ALL="ALL"

def callback_for_args(args):
    if args.fdr_filter is not None:
        filter = PB8KFileInfo.FDRFilter(args.fdr_filter)
    else:
        filter = None

    if args.model_weight_type == BETA:
        if not args.variance_path or not args.sample_size:
            raise RuntimeError("Insufficient parameters for beta model")
        row_from_comps = PB8KFileInfo.BetaRowFromComps(args.variance_path, args.sample_size)
    elif args.model_weight_type == ZSCORE:
        row_from_comps = PB8KFileInfo.ZScoreRowFromComps()
    else:
        raise RuntimeError("Wrong building type %s", args.model_building_type)

    if args.model_building_criteria == ALL:
        process_rows = PB8KFileInfo.process_all_rows
    elif args.model_building_criteria == BEST:
        process_rows = PB8KFileInfo.keep_best_row
    else:
        raise RuntimeError("Wrong building criteria: %s", args.model_building_criteria)

    callback = PB8KFileInfo.PB8KFileCallback(row_from_comps, process_rows, filter)
    return callback


def run(args):
    if os.path.exists(args.model_output_path):
        logging.info("DB already there, delete it if you want it done again")
        return

    folder = os.path.split(args.model_output_path)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    callback = callback_for_args(args)
    PB8KFileInfo.parse_input_file(args.model_output_path, args.pb8k_input_path, args.gencode_input_path, callback)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Build Intron Quantification and Geuvadis phenotype to transcriptome model pipeline input")

    parser.add_argument("--pb8k_input_path",
                        help="path to PB8K data files",
                        default="data/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt.gz")

    parser.add_argument("--gencode_input_path",
                        help="path to gencode file",
                        default="data/gencode.v22.annotation.gtf.gz")

    parser.add_argument("--model_output_path",
                        help="path to model file",
                        default="results/PRS_PB8K_ZSCORE_ALL.db")

    parser.add_argument("--variance_path",
                        help="Optional: snp`variance file",
                        default=None)
#                        default="data/VAR_TGF_EUR_PB8K.txt.gz")

    parser.add_argument("--fdr_filter",
                        help="Optional false discovery ratio filter",
                        type=float,
                        default=None)

    parser.add_argument("--sample_size",
                        help="data sample size",
                        type=int,
                        default=5311)

    parser.add_argument("--model_building_criteria",
                        help="Type of PRS model. [BEST/ALL] snps in a gene.",
                        default=ALL)

    parser.add_argument("--model_weight_type",
                        help="Type of PRS model wright. [ZSCORE/BETA]; BETA requires variance and sample size parameter",
                        default=ZSCORE)

    parser.add_argument("--expect_throw",
                    help="Debug mode for incomplete gzipped file",
                    action="store_true",
                    default=False)


    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    run(args)