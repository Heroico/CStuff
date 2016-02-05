#! /usr/bin/env python
import logging
import gzip
import Logging
import WeightDBUtilities

def parse_file(path, callback):
    with gzip.open(path) as file:
        try:
            for i, line in enumerate(file):
                callback(i, line.strip().split())
        except IOError:
            logging.info("Expected bad file, and that's what we got: %s" % (path, ))
        except Exception as e:
            raise e

#snp dict
class SNPTF(object):
    Chr=0
    Pos=1
    VariantID=2
    Ref_b37=3
    Alt=4
    RS_ID_dbSNP135_original_VCF=5
    RS_ID_dbSNP142_CHG37p13=6
    Num_alt_per_site=7


class BuildSnpDictCallback(object):
    def __init__(self, weight_db):
        self.weight_db = weight_db
        self.snps = {}

    def __call__(self, i, comps):
        if i==0:
            return

        if len(comps) != 8:
            return

        snp = comps[SNPTF.RS_ID_dbSNP142_CHG37p13]
        if snp in self.weight_db.genes_for_an_rsid:
            self.snps[snp] = (comps[SNPTF.VariantID], comps[SNPTF.Ref_b37], comps[SNPTF.Alt])


def build_snp_dict(path, weight_db_logic):
    build_snp = BuildSnpDictCallback(weight_db_logic)
    parse_file(path, build_snp)
    return build_snp.snps

# genotype
class BuildGeneExpressionCallback(object):
    def __init__(self, weight_db, snp_dict):
        self.weight_db = weight_db
        self.snp_dict = snp_dict
        self.collected = {}

    def __call__(self, i, comps):
        if i == 0:
            comps = comps[3:] #individuals
            return
        #snp geno7

def build_gene_expression(path, weight_db, snp_dict):
    callback = BuildGeneExpressionCallback(weight_db, snp_dict)
    parse_file(path, callback)
    return callback.collected

class BuildExpressionWork(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        logging.info("Loading %s", self.args.weight_db)
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.args.weight_db)

        logging.info("Building snp dict from %s", self.args.gtex_snp)
        snp_dict = build_snp_dict(self.args.gtex_snp, weight_db_logic)

        logging.info("Building gene expression")
        gene_expression = build_gene_expression(self.args.gtex_geno, weight_db_logic, snp_dict)

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