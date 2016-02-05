import Utilities

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
    Utilities.parse_file(path, build_snp)
    return build_snp.snps
