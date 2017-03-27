#! /usr/bin/env python
import gzip

SUPPL_FOLDER = "/scratch/abarbeira3/test/data"
EQTL_FOLDER= "/group/im-lab/nas40t2/Data/dbGaP/GTEx/V6p-gtexportal/GTEx_Analysis_v6p_eQTL/"

#GENCODE="/group/im-lab/nas40t2/abarbeira/data/GTEx/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt.gz"
GENCODE="/home/numa/Documents/Projects/data/GTEx/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt.gz"

def load_gencode(path):
    results={}
    with gzip.open(path) as file:
        for line in file:
            if "#" in line: continue
            comps = line.strip().split()
            variant=comps[2]
            rsid=comps[6]
            results[variant] = rsid
    return results

results = load_gencode(GENCODE)
from IPython import embed; embed()


