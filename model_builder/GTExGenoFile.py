import logging
import Utilities
import GTExSNPFile
import re

# Gtex geno file format
class GPF(object):
    SNP=0
    GENE=1
    BETA=2
    T_STAT=3
    P_VALUE=4

# genotype
class BuildGeneExpressionCallback(object):
    def __init__(self, weight_db, snp_dict):
        self.weight_db = weight_db
        self.snp_dict = snp_dict
        self.collected = {}
        self.individuals = []

    def __call__(self, i, comps):
        if i == 0:
            regexp = re.compile('^[\w-]+$')
            self.individuals = [x for x in comps[1:] if regexp.match(x) is not None and not "Id" in x] #individuals
            logging.info("Individuals: %d", len(self.individuals))
            return
        #snp geno7
        variant = comps[0]
        if not variant in self.snp_dict:
            logging.log(5, "Missing variant %s in snp_dict", variant)
            return

        variant_comps = variant.split("_")
        ref_allele = variant_comps[2]
        eff_allele = variant_comps[3]

        snp_info = self.snp_dict[variant]
        rsid = snp_info[GTExSNPFile.SNP_E_F.RSID]
        if not rsid in self.weight_db.genes_for_an_rsid:
            logging.log(6, "Missing variant %s, rsid %s in weight_db", variant, rsid)
            return

        dosages = comps[1:]
        if len(dosages) != len(self.individuals):
            logging.log(5,"Bad dosage line: %s", str(comps))
            return

        for gene in self.weight_db.genes_for_an_rsid[rsid]:
            weight = self.weight_db.weights_by_gene_name[gene][rsid]
            if {ref_allele, eff_allele} != {weight.ref_allele, weight.eff_allele}:
                logging.log(6, "Wrong alleles for variant %s, rsid %s in weight_db: %s,%s  %s, %s",
                            variant, rsid, ref_allele, eff_allele, weight.ref_allele, weight.eff_allele)
                continue

            if not gene in self.collected:
                self.collected[gene] = {x:0 for x in self.individuals}

            flip = ref_allele == weight.eff_allele and eff_allele == weight.ref_allele

            exp = self.collected[gene]

            for i, dosage in enumerate(dosages):
                d = float(dosage)
                exp[self.individuals[i]] += (2-d if flip else d) * weight.weight


def build_gene_expression(path, weight_db, snp_dict):
    callback = BuildGeneExpressionCallback(weight_db, snp_dict)
    Utilities.parse_file(path, callback)
    return callback.collected, callback.individuals