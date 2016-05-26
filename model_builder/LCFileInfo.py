import Utilities
import GTExSNPFile

# Using tuples instead of objects for memory consumption reasons!

#At the actual file, snp is snp[chr]_[position]
class LCTF(object):
    INTRON_CLUSTER=0
    SNP=1
    PVALUE=2
    BETA=3
    STE=4
    REFERENCE_ALLELE=5
    ALTERNATE_ALLELE=6

class LCTFBasicLoadCallback(object):
    def __init__(self, format=LCTF):
        self.format = format
        self.results = []

    def __call__(self, i, comps):
        f = self.format
        e = (comps[f.INTRON_CLUSTER], comps[f.SNP], comps[f.PVALUE], comps[f.BETA], comps[f.STE], comps[f.REFERENCE_ALLELE],comps[f.ALTERNATE_ALLELE],)
        self.results.append(comps)

def load_lctf(file_path, format=LCTF):
    c = LCTFBasicLoadCallback(format)
    Utilities.parse_file(file_path, c)
    return c.results

def fill_snp_id_from_gtex(entries, gtex_file_path, lctf_format=LCTF, gtex_format=GTExSNPFile.SNPTF, expect_throw=False):
    class TheCallback(object):
        def __init__(self, entries, lctf_format, gtex_format):
            self.entries = {e[lctf_format.SNP]:e for e in entries}
            self.lctf_format = lctf_format
            self.gtex_format = gtex_format
            self.results = []

        def __call__(self, i, line):
            if i==0:
                return

            SNPTF = self.gtex_format
            VARIANT_ID_F = GTExSNPFile.VARIANT_ID_F
            rf = self.lctf_format

            variant = line[SNPTF.VariantID]
            variant_comps = variant.split("_")
            ref_a = variant_comps[VARIANT_ID_F.REF_ALLELE]
            alt_a = variant_comps[VARIANT_ID_F.ALT_ALLELE]
            pos = line[SNPTF.Pos]

            #we chose lctf snp_[chr]_pos as key
            snp_pos = "_".join(["snp", variant_comps[VARIANT_ID_F.CHR],pos])

            # if we cannot get an rsid, we'lll los e it
            if not snp_pos in self.entries:
                return

            entry = self.entries[snp_pos]
            if {ref_a, alt_a} != {entry[rf.REFERENCE_ALLELE], entry[rf.ALTERNATE_ALLELE]}:
                return

            rsid = line[SNPTF.RS_ID_dbSNP142_CHG37p13]
            new_entry = (entry[rf.INTRON_CLUSTER], rsid, entry[rf.PVALUE], entry[rf.BETA], entry[rf.STE], entry[rf.REFERENCE_ALLELE], entry[rf.ALTERNATE_ALLELE],)
            self.results.append(new_entry)
            del self.entries[snp_pos]

    c = TheCallback(entries, lctf_format, gtex_format)
    Utilities.parse_file(gtex_file_path, c, expect_throw=expect_throw)
    return c.results

def entry_to_weight_db(entries, name_from_entry, weight_from_entry, f=LCTF, wf=Utilities.WDBIF):
    db_entries = []
    for e in entries:
        name = name_from_entry(e, f)
        w = weight_from_entry(e, f)
        # I think alleles should be inverted
        db_entry = (e[f.SNP], name, name, w, e[f.REFERENCE_ALLELE], e[f.ALTERNATE_ALLELE])
        db_entries.append(db_entry)

    r = {}
    for db_e in db_entries:
        gene = db_e[wf.GENE]
        if not gene in r:
            r[gene] = []
        rows = r[gene]
        rows.append(db_e)

    return r

def cluster_name_from_entry(entry, f):
    cluster = entry[f.INTRON_CLUSTER]
    name = "clu_"+cluster.split("clu_")[1]
    return name

def zscore_from_entry(entry, f):
    beta = float(entry[f.BETA])
    se = float(entry[f.STE])
    w = str(beta/se)
    return w

def beta_from_entry(entry, f):
    return entry[f.BETA]

def entries_to_zscore_weight_db_by_cluster(entries, f=LCTF, wf=Utilities.WDBIF):
    r = entry_to_weight_db(entries, cluster_name_from_entry, zscore_from_entry, f, wf)
    return r

def entries_to_beta_weight_db_by_cluster(entries, f=LCTF, wf=Utilities.WDBIF):
    r = entry_to_weight_db(entries, cluster_name_from_entry, beta_from_entry, f, wf)
    return r