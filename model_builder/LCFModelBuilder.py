#! /usr/bin/env python
import logging
import Logging
import LCFileInfo
import Utilities

#ran with:
# nohup python LCFModelBuilder.py --lctf_path /group/im-lab/nas40t2/Data/SummaryResults/protected/Leafcutter/sQTLs_geuv_summary.gen.txt.gz --gtex_snp_path /group/im-lab/nas40t2/Data/dbGaP/GTEx/V6-gtexportal/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz --output_file /group/im-lab/nas40t2/abarbeira/Projects/model_dbs/PRS_sQTL_geuv.db > lctf.log 2>&1 &

NAME_METHODS = {
    "CLUSTER":LCFileInfo.cluster_name_from_entry,
    "INTRON":LCFileInfo.intron_name_from_entry
}

WEIGHT_METHODS = {
    "BETA":LCFileInfo.beta_from_entry,
    "ZSCORE":LCFileInfo.zscore_from_entry
}

def process(args):
    logging.info("Loading LCTF %s", args.lctf_path)
    lctf_entries = LCFileInfo.load_lctf(args.lctf_path)

    logging.info("Filling snp rsid from %s", args.gtex_snp_path)
    lctf_entries = LCFileInfo.fill_snp_id_from_gtex(lctf_entries, args.gtex_snp_path, expect_throw=args.expect_throw)

    logging.info("Building db entries")
    #This will drop the previous stuff.
    weight_method= WEIGHT_METHODS[args.weight_type]
    name_method = NAME_METHODS[args.functional_unit_type]
    lctf_entries = LCFileInfo.entries_to_weight_db(lctf_entries, name_method, weight_method)

    logging.info("Building db")
    connection = Utilities.connect(args.output_file)

    Utilities.setup_db(connection)

    logging.info("Inserting entries")
    Utilities.insert_entries(connection, lctf_entries)

    logging.info("Ran successfully")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Build weight db model")

    parser.add_argument("--lctf_path",
                        help="path to lctf file",
                        default="data/sQTLs_geuv_summary.gen.txt.gz")

    parser.add_argument("--gtex_snp_path",
                        help="path to GTEx snp file",
                        default="data/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz")

    parser.add_argument("--output_file",
                        help="output where db will be saved",
                        default="results/PRS_sQTLs_geuv.db")

    parser.add_argument("--functional_unit_type",
                        help="INTRON or CLUSTER",
                        default="INTRON")

    parser.add_argument("--weight_type",
                        help="ZSCORE or BETA",
                        default="BETA")

    parser.add_argument("--expect_throw",
                    help="Debug mode for incomplete gzipped file",
                    action="store_true",
                    default=False)


    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    process(args)