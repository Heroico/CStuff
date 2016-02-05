import Utilities
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
    Utilitites.parse_file(path, callback)
    return callback.collected