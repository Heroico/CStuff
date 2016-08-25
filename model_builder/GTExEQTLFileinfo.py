import logging
import os
import math
import numpy
import qvalue
import Utilities
import GencodeFile
import GTExSNPFile

class GTEXEQTLF(object):
    VARIANT_ID = 0
    GENE_ID = 1
    TSS_DISTANCE = 2
    PVAL_NOMINAL = 3
    SLOPE = 4
    SLOPE_SE = 5
    SLOPE_FPKM = 6
    SLOPE_FPKM_SE = 7
    PVAL_NOMINAL_THRESHOLD = 8
    MIN_PVAL_NOMINAL = 9
    PVAL_BETA = 10

    HEADER="variant_id\tgene_id\ttss_distance\tpval_nominal\tslope\tslope_se\tslope_fpkm\tslope_fpkm_se\tpval_nominal_threshold\tmin_pval_nominal\tpval_beta"

class ZScoreRowFromComps(object):
    def __init__(self, tf=GTEXEQTLF):
        self.TF = tf

    def __call__(self, i, comps):
        gene_id = comps[self.TF.GENE_ID]
        snp = comps[self.TF.VARIANT_ID]
        beta = comps[self.TF.SLOPE]
        beta_se = comps[self.TF.SLOPE_SE]
        zscore = str(float(beta)/float(beta_se))
        p = comps[self.TF.PVAL_NOMINAL]
        row = (snp, gene_id, None, None, None, zscore, None, None, p, None)
        return row

class BetaRowFromComps(object):
        def __init__(self, tf=GTEXEQTLF):
            self.TF = tf

        def __call__(self, i, comps):
            gene_id = comps[self.TF.GENE_ID]
            snp = comps[self.TF.VARIANT_ID]
            beta = comps[self.TF.SLOPE]
            p = comps[self.TF.PVAL_NOMINAL]
            row = (snp, gene_id, None, None, None, beta, None, None, p, None)
            return row

class GTExEQTLFileCallback(object):
    def __init__(self, row_from_comps, process_row, TF=GTEXEQTLF):
        self.genes = {}
        self.row_from_comps = row_from_comps
        self.process_row = process_row
        self.TF = TF
        self.pvalues = []

    def __call__(self, i, comps):
        row = self.row_from_comps(i, comps)
        self.process_row(row, self.genes, Utilities.WDBIF.GENE)
        self.pvalues.append(float(row[Utilities.WDBIF.GENE_PVALUE]))

def fix_rows(genes, pvalues):
    pvalues = numpy.array(pvalues)
    _qvalues = qvalue.estimate(pvalues)
    _qvalues = {pvalues[i]:q for i,q in enumerate(_qvalues)}
    F = Utilities.WDBIF
    keys = genes.keys()
    for gene in keys:
        rows = genes[gene]
        input = [(math.fabs(float(r[F.WEIGHT])), float(r[F.GENE_PVALUE]), _qvalues[float(r[F.GENE_PVALUE])], r[F.SNP]) for k,r in rows.iteritems()]
        w = sum(map(lambda x: x[0], input))
        if w == 0:
            logging.info("Cannot handle null weight for %s", gene)
            a_p = None
            a_q = None
        else:
            a_p = str(sum(map(lambda x: math.fabs(x[0])*x[1], input)) / w)
            a_q = str(sum(map(lambda x: math.fabs(x[0])*x[2], input)) / w)
        n = str(len(input))
        genes[gene] = [(r[F.SNP], r[F.GENE], r[F.GENE_NAME], r[F.REFERENCE_ALLELE], r[F.EFFECT_ALLELE], r[F.WEIGHT], n, r[F.GENE_R2], a_p, a_q) for k,r in rows.iteritems()]
    return genes

def parse_input_file(db_output_path, gtex_pheno_input_path, gencode_input_path, gtex_snp_path, gtex_callback):
    logging.info("Opening GTEx pheno file %s", os.path.basename(os.path.normpath(gtex_pheno_input_path)))
    file_iterator = Utilities.CSVFileIterator(gtex_pheno_input_path, delimiter="\t", header=GTEXEQTLF.HEADER, compressed=True)
    file_iterator.iterate(gtex_callback)
    logging.info("%d found at GTEx file", len(gtex_callback.genes))

    logging.info("Opening gencode file")
    def gencode_fixed_row(gencode, row):
        F = Utilities.WDBIF
        return (row[F.SNP], row[F.GENE], gencode.gene_name, row[F.REFERENCE_ALLELE], row[F.EFFECT_ALLELE], row[F.WEIGHT], row[F.N_SNP], row[F.GENE_R2], row[F.GENE_PVALUE], row[F.GENE_QVALUE])

    class FixGeneCallback(object):
        def __init__(self, genes):
            self.genes = genes
            self.selected = {}

        def __call__(self, gencode):
            if gencode.gene_id in self.genes:
                rows = self.genes[gencode.gene_id]
                self.selected[gencode.gene_id] = {k:gencode_fixed_row(gencode, row) for k,row in rows.iteritems()}

    gencode_callback = FixGeneCallback(gtex_callback.genes)
    GencodeFile.parse_gencode_file(gencode_input_path, gencode_callback)
    genes = gencode_callback.selected
    logging.info("%d survived after gencode file", len(genes))
    pvalues = gtex_callback.pvalues
    del gencode_callback
    del gtex_callback

    logging.info("Opening GTEX Snp file")
    def snp_fixed_row(row, rsid, ref_allele, eff_allele):
        F = Utilities.WDBIF
        return (rsid, row[F.GENE], row[F.GENE_NAME], ref_allele, eff_allele, row[F.WEIGHT], row[F.N_SNP], row[F.GENE_R2], row[F.GENE_PVALUE], row[F.GENE_QVALUE])

    class FixSNPCallback(object):
        def __init__(self, genes):
            self.genes = genes
            self.selected = {}
            self.variant_to_gene = {}
            for gene, rows in self.genes.iteritems():
                for variant, row in rows.iteritems():
                    self.variant_to_gene[variant] = gene

        def __call__(self, i, comps):
            if i == 0:
                return
            F = GTExSNPFile.SNPTF
            variant = comps[F.VariantID]
            ref = comps[F.Ref_b37]
            eff = comps[F.Alt]
            snp = comps[F.RS_ID_dbSNP142_CHG37p13]
            if variant in self.variant_to_gene:
                gene = self.variant_to_gene[variant]
                row = self.genes[gene][variant]
                if not gene in self.selected:
                    self.selected[gene] = {}
                self.selected[gene][snp] = snp_fixed_row(row, snp, ref, eff)

    snp_callback = FixSNPCallback(genes)
    snp_iterator = Utilities.CSVFileIterator(gtex_snp_path, delimiter="\t", compressed=True) #header in the file is just wrong., header=GTExSNPFile.SNPTF.HEADER, compressed=True)
    snp_iterator.iterate(snp_callback)
    genes = snp_callback.selected
    del snp_callback.selected

    logging.info("Fixing rows")
    fix_rows(genes, pvalues)

    logging.info("Saving entries")
    connection = Utilities.connect(db_output_path)
    Utilities.setup_db(connection)
    Utilities.insert_entries(connection, genes)
