#! /usr/bin/env python
import logging
import os

import Logging
import GTExEQTLFileinfo
import Utilities

BETA="BETA"
ZSCORE="ZSCORE"


BEST="BEST"
ALL="ALL"

def callback_for_args(args):
    # if args.fdr_filter is not None:
    #     filter = PB8KFileInfo.FDRFilter(args.fdr_filter)
    # else:
    #     filter = None
    MODEL_TYPE = {
        BETA:GTExEQTLFileinfo.BetaRowFromComps(),
        ZSCORE:GTExEQTLFileinfo.ZScoreRowFromComps()
    }
    if args.model_weight_type not in MODEL_TYPE:
        raise RuntimeError("Wrong model type: %s", args.model_weight_type)
    row_from_comps = MODEL_TYPE[args.model_weight_type]

    MODEL_CRITERIA = {
        ALL:Utilities.process_all_rows,
        BEST:Utilities.keep_best_row
    }
    if args.model_building_criteria not in MODEL_CRITERIA:
        raise RuntimeError("Wrong building criteria: %s", args.model_building_criteria)
    process_row = MODEL_CRITERIA[args.model_building_criteria]

    callback = GTExEQTLFileinfo.GTExEQTLFileCallback(row_from_comps, process_row, filter)
    return callback

def run(args):
    if os.path.exists(args.model_output_path):
        logging.info("DB already there, delete it if you want it done again")
        return

    folder = os.path.split(args.model_output_path)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    the_callback = callback_for_args(args)
    GTExEQTLFileinfo.parse_input_file(db_output_path=args.model_output_path,
                                      gtex_pheno_input_path=args.gtex_eqtl_input_path,
                                      gtex_snp_path=args.gtex_snp,
                                      gencode_input_path=args.gencode_input_path,
                                      gtex_callback=the_callback)
    logging.info("Ran successfully")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Build Polygenic risk score model")

    parser.add_argument("--gtex_eqtl_input_path",
                        help="path to PB8K data files",
                        default="data/GTEx_Analysis_v6p_eQTL/Whole_Blood_Analysis.v6p.signif_snpgene_pairs.txt.gz")

    parser.add_argument("--gencode_input_path",
                        help="path to gencode file",
                        default="data/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz")

    parser.add_argument("--gtex_snp",
                        help= "snp metadata file",
                        default="data/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz")

    parser.add_argument("--model_output_path",
                        help="path to model file",
                        default="results/Whole_Blood_Beta_BEST.db")

    parser.add_argument("--fdr_filter",
                        help="Optional false discovery ratio filter",
                        type=float,
                        default=None)

    parser.add_argument("--model_building_criteria",
                        help="Type of PRS model. [BEST/ALL] snps in a gene.",
                        default=BEST)

    parser.add_argument("--model_weight_type",
                        help="Type of PRS model wright. [ZSCORE/BETA]; BETA requires variance and sample size parameter",
                        default=BETA)

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