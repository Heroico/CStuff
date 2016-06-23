#! /usr/bin/env python
import logging
import gzip
import os
import Logging
import WeightDBUtilities
import GTExGenoFile
import GTExSNPFile

def save_expression(path, gene_expression, individuals):
    sorted_individuals = sorted(individuals)
    with gzip.open(path, "wb") as file:
        header = "\t".join( ["Id"]+sorted_individuals ) + "\n"
        file.write(header)

        keys = sorted(gene_expression.keys())
        for key in keys:
            entries = [str(gene_expression[key][individual]) for individual in sorted_individuals]
            line = "\t".join( [key] + entries) + "\n"
            file.write(line)


class BuildExpressionWork(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        if os.path.exists(self.args.output):
            logging.info("%s already exists. Delete it if you want it done again", self.args.output)
            return

        logging.info("Loading %s", self.args.weight_db)
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.args.weight_db)

        logging.info("Building snp dict from %s", self.args.gtex_snp)
        snp_dict = GTExSNPFile.build_snp_dict(self.args.gtex_snp, weight_db_logic)

        logging.info("Building gene expression")
        gene_expression, individuals = GTExGenoFile.build_gene_expression(self.args.gtex_geno, weight_db_logic, snp_dict)

        logging.info("Saving gene expression")
        save_expression(self.args.output, gene_expression, individuals)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Build gene expresion from weight db")

    parser.add_argument("--weight_db",
                        help="path to weight db",
                        default="results/PRS_PB8K_BETA_ALL.db")

    parser.add_argument("--gtex_geno",
                        help="Path to gtex genotype",
                        default="data/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices.tar.gz")

    parser.add_argument("--gtex_snp",
                        help= "snp metadata file",
                        default="data/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz")

    parser.add_argument("--output",
                        help="where to save expression",
                        default="results/PRS_PB8K_BETA_ALL_expression.txt.gz")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()
    Logging.configureLogging(int(args.verbosity))

    work = BuildExpressionWork(args)
    work.run()