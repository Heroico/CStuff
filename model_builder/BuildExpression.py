#! /usr/bin/env python
import logging
import Logging
import WeightDBUtilities
import GTExGenoFile
import GTExSNPFile

class BuildExpressionWork(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        logging.info("Loading %s", self.args.weight_db)
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.args.weight_db)

        logging.info("Building snp dict from %s", self.args.gtex_snp)
        snp_dict = GTExSNPFile.build_snp_dict(self.args.gtex_snp, weight_db_logic)

        logging.info("Building gene expression")
        gene_expression = GTExGenoFile.build_gene_expression(self.args.gtex_geno, weight_db_logic, snp_dict)

if __name__ == "__main__":
    class Args(object):
        def __init__(self):
            self.weight_db = "results/PB8K.db"
            self.gtex_geno = "data/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices.tar.gz"
            self.gtex_snp = "data/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz"
            self.verbosity = 10

    args = Args()
    Logging.configureLogging(int(args.verbosity))

    work = BuildExpressionWork(args)
    work.run()