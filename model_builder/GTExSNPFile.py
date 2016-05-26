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

class SNP_E_F(object):
    RSID=0
    VARIANT=1
    REF_ALLELE=0
    EFF_ALLELE=1

class VARIANT_ID_F(object):
    CHR=0
    POS=1
    REF_ALLELE=2
    ALT_ALLELE=3
    VERSION=4

class BuildSnpDictCallback(object):
    def __init__(self, weight_db):
        self.weight_db = weight_db
        self.snps = {}

    def __call__(self, i, comps):
        if i==0:
            return

        if len(comps) != 8:
            return

        rsid = comps[SNPTF.RS_ID_dbSNP142_CHG37p13]
        variant = comps[SNPTF.VariantID]
        if rsid in self.weight_db.genes_for_an_rsid:
            self.snps[variant] = (rsid, variant, comps[SNPTF.Ref_b37], comps[SNPTF.Alt])


def build_snp_dict(path, weight_db_logic):
    build_snp = BuildSnpDictCallback(weight_db_logic)
    Utilities.parse_file(path, build_snp)
    return build_snp.snps
