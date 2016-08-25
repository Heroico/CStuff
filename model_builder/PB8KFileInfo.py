import logging
import math
import Utilities
import VarianceFile
import GencodeFile

class TF1(object):
    PValue=0
    SNPName=1
    SNPChr=2
    SNPChrPos=3
    ProbeName=4
    ProbeChr=5
    ProbeCenterChrPos=6
    CisTrans=7
    SNPType=8
    AlleleAssessed=9
    OverallZScore=10
    DatasetsWhereSNPProbePairIsAvailableAndPassesQC=11
    DatasetsZScores=12
    HUGO=13
    FDR=14

    HEADER="PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tHUGO\tFDR"

def alleles(TF, comps):
    alleles = {a for a in comps[TF.SNPType].split("/")}

    #TODO: figure this out
    assessed = comps[TF.AlleleAssessed]
    for a in alleles:
        if a == assessed:
            effect_allele = a
        else:
            reference_allele = a
    return reference_allele, effect_allele

class FDRFilter(object):
    def __init__(self, fdr_filter, TF=TF1):
        self.TF = TF
        self.fdr_filter = fdr_filter

    def __call__(self, comps):
        fdr = comps[self.TF.FDR]
        if float(fdr) > self.fdr_filter:
            snp = comps[self.TF.SNPName]
            logging.log(9,"Snp %s doesn't pass fdr filter: %s", snp, fdr)
            return False
        return True

# Function-like object to act as row builder
# Zscore as weight, tentatively loading Pvalue and FDR even though it might discarded
# See Utilities.WDBIF format for output
class ZScoreRowFromComps(object):
    def __init__(self, tf=TF1):
        self.TF = tf

    def __call__(self, i, comps, gene_name):
        snp = comps[self.TF.SNPName]
        reference_allele, effect_allele = alleles(self.TF, comps)
        zscore = comps[self.TF.OverallZScore]
        p = comps[self.TF.PValue]
        q = comps[self.TF.FDR]
        row = (snp, None, gene_name, reference_allele, effect_allele, zscore, "NA", "NA", p, q)
        return row

# Function-like object to act as row builder
# BETA as weight, tentatively loading Pvalue and FDR even though it might discarded
# See Utilities.WDBIF format for output
class BetaRowFromComps(object):
    def __init__(self, variance_file_path, sample_size, TF=TF1):
        self.TF = TF
        logging.info("Opening variance file")
        self.vars = VarianceFile.load_variance(variance_file_path)
        self.sample_size = sample_size

    def __call__(self, i, comps, gene_name):
        snp = comps[self.TF.SNPName]
        reference_allele, effect_allele = alleles(self.TF, comps)
        zscore = comps[self.TF.OverallZScore]
        if not snp in self.vars:
            return None
        v = self.vars[snp]
        std = math.sqrt(v/self.sample_size)
        b = str(float(zscore)*std)
        p = comps[self.TF.PValue]
        q = comps[self.TF.FDR]
        row = (snp, None, gene_name, reference_allele, effect_allele, b, "NA", "NA", p, q)
        return row

# Master callback for each row in PB8K file
class PB8KFileCallback(object):
    def __init__(self, row_from_comps, process_row, row_filter, TF=TF1):
        self.genes = {}
        self.row_from_comps = row_from_comps
        self.process_row = process_row
        self.row_filter = row_filter
        self.TF = TF

    def __call__(self, i, comps):
        if self.row_filter:
            passes_filter = self.row_filter(comps)
            if not passes_filter:
                return
        gene_name = comps[self.TF.HUGO]
        if "," in gene_name:
            multiple_genes = gene_name.split(",")
            for g in multiple_genes:
                row = self.row_from_comps(i, comps, g)
                if row:
                    self.process_row(row, self.genes, Utilities.WDBIF.GENE_NAME)
        else:
            row = self.row_from_comps(i, comps, gene_name)
            if row:
                self.process_row(row, self.genes, Utilities.WDBIF.GENE_NAME)

def fix_row(genes):
    F = Utilities.WDBIF
    keys = genes.keys()
    for gene in keys:
        rows = genes[gene]
        input = [(math.fabs(float(r[F.WEIGHT])), float(r[F.GENE_PVALUE]), float(r[F.GENE_QVALUE]), r[F.SNP]) for r in rows]
        w = sum(map(lambda x: x[0], input))
        if w == 0:
            logging.info("Cannot handle null variance for %s", gene)
            a_p = None
            a_q = None
        else:
            a_p = str(sum(map(lambda x: math.fabs(x[0])*x[1], input)) / w)
            a_q = str(sum(map(lambda x: math.fabs(x[0])*x[2], input)) / w)
        n = str(len(input))
        genes[gene] = [(r[F.SNP], r[F.GENE], r[F.GENE_NAME], r[F.REFERENCE_ALLELE], r[F.EFFECT_ALLELE], r[F.WEIGHT], n, r[F.GENE_R2], a_p, a_q) for r in rows]
    return genes

def parse_input_file(db_output_path, pheno_input_path, gencode_input_path,  pb8k_callback):
    logging.info("Opening PB8K pheno file")
    file_iterator = Utilities.CSVFileIterator(pheno_input_path, delimiter="\t", header=TF1.HEADER, compressed=True)
    file_iterator.iterate(pb8k_callback)

    logging.info("Opening gencode file")
    def fixed_row(gencode, row):
        F = Utilities.WDBIF
        return (row[F.SNP], gencode.ensemble_version, row[F.GENE_NAME], row[F.REFERENCE_ALLELE], row[F.EFFECT_ALLELE], row[F.WEIGHT], row[F.N_SNP], row[F.GENE_R2], row[F.GENE_PVALUE], row[F.GENE_QVALUE])

    class FixGeneCallback(object):
        def __init__(self, genes):
            self.genes = genes
            self.selected = {}

        def __call__(self, gencode):
            if gencode.name in self.genes:
                rows = self.genes[gencode.name]
                self.selected[gencode.name] = [fixed_row(gencode, row) for k,row in rows.iteritems()]

    gencode_callback = FixGeneCallback(pb8k_callback.genes)
    GencodeFile.parse_gencode_file(gencode_input_path, gencode_callback)
    genes = gencode_callback.selected

    logging.info("Fixing rows")
    genes = fix_row(genes)

    logging.info("Saving entries")
    connection = Utilities.connect(db_output_path)
    Utilities.setup_db(connection)
    Utilities.insert_entries(connection, genes)
