import scipy.stats as stats
import logging
import math
import os
import numpy

import VarianceFile
import Utilities
import GTExSNPFile
import qvalue

class FTF(object):
    GENE_SYMBOL=0
    ENSEMBL_ID=1
    INDEX_E_SNP=2
    P_VALUE=3
    ADJUSTED_P_VALUE=4

    HEADER="Gene symbol,ENSEMBL ID,Index eSNP,p-value,adjusted p-value"

# Function-like object to act as row builder
# Zscore as weight, tentatively loading Pvalue and FDR even though it might discarded
# See Utilities.WDBIF format for output
class ZScoreRowFromComps(object):
    def __init__(self, tf=FTF):
        self.TF = tf

    def __call__(self, i, comps):
        gene_name = comps[self.TF.GENE_SYMBOL]
        gene_id = comps[self.TF.ENSEMBL_ID]
        snp = comps[self.TF.INDEX_E_SNP]
        p = comps[self.TF.P_VALUE]
        zscore = p_to_z(p)
        row = (snp, gene_id, gene_name, "NA", "NA", zscore, "NA", "NA", p, "NA")
        return row

def p_to_z(p):
    zscore = -stats.norm.ppf(float(p)/2)
    return zscore

# Function-like object to act as row builder
# BETA as weight, tentatively loading Pvalue and FDR even though it might discarded
# See Utilities.WDBIF format for output
class BetaRowFromComps(object):
    def __init__(self, variance_file_path, sample_size, TF=FTF):
        self.TF = TF
        logging.info("Opening variance file")
        self.vars = VarianceFile.load_variance(variance_file_path)
        self.sample_size = sample_size

    def __call__(self, i, comps):
        gene_name = comps[self.TF.GENE_SYMBOL]
        gene_id = comps[self.TF.ENSEMBL_ID]
        snp = comps[self.TF.INDEX_E_SNP]
        p = comps[self.TF.P_VALUE]
        zscore = p_to_z(p)
        if not snp in self.vars:
            return None
        v = self.vars[snp]
        std = math.sqrt(v/self.sample_size)
        b = str(float(zscore)*std)
        row = (snp, gene_id, gene_name, "NA", "NA", b, "NA", "NA", p, "NA")
        return row

class FranzenEQTLFileCallback(object):
    def __init__(self, row_from_comps, process_row, TF=FTF):
        self.genes = {}
        self.row_from_comps = row_from_comps
        self.process_row = process_row
        self.TF = TF
        self.pvalues = []

    def __call__(self, i, comps):
        if i <2:
            return
        row = self.row_from_comps(i, comps)
        if row is not None:
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

def parse_input_file(db_output_path, franzen_pheno_input_path, gtex_snp_path, franzen_callback):
        logging.info("Opening Frazen pheno file %s", os.path.basename(os.path.normpath(franzen_pheno_input_path)))
        file_iterator = Utilities.CSVFileIterator(franzen_pheno_input_path, delimiter=",", header=FTF.HEADER, compressed=True)
        file_iterator.iterate(franzen_callback)
        logging.info("%d found at Franzen file", len(franzen_callback.genes))
        selected = franzen_callback.genes
        pvalues = franzen_callback.pvalues
        del franzen_callback

        # logging.info("Opening gencode file")
        #
        # def gencode_fixed_row(gencode, row):
        #     F = Utilities.WDBIF
        #     return (
        #     row[F.SNP], row[F.GENE], gencode.gene_name, row[F.REFERENCE_ALLELE], row[F.EFFECT_ALLELE], row[F.WEIGHT],
        #     row[F.N_SNP], row[F.GENE_R2], row[F.GENE_PVALUE], row[F.GENE_QVALUE])
        #
        # class FixGeneCallback(object):
        #     def __init__(self, genes):
        #         self.genes = genes
        #         self.selected = {}
        #
        #     def __call__(self, gencode):
        #         if gencode.gene_id in self.genes:
        #             rows = self.genes[gencode.gene_id]
        #             self.selected[gencode.gene_id] = {k: gencode_fixed_row(gencode, row) for k, row in rows.iteritems()}
        #
        # gencode_callback = FixGeneCallback(gtex_callback.genes)
        # GencodeFile.parse_gencode_file(gencode_input_path, gencode_callback)
        # genes = gencode_callback.selected
        # logging.info("%d survived after gencode file", len(genes))
        # pvalues = gtex_callback.pvalues
        # del gencode_callback
        # del gtex_callback
        #
        logging.info("Opening GTEX Snp file")
        def snp_fixed_row(row, rsid, ref_allele, eff_allele):
            F = Utilities.WDBIF
            return ( rsid, row[F.GENE], row[F.GENE_NAME], ref_allele, eff_allele, row[F.WEIGHT], row[F.N_SNP], row[F.GENE_R2], row[F.GENE_PVALUE], row[F.GENE_QVALUE])

        class FixSNPCallback(object):
            def __init__(self, genes):
                self.genes = genes
                self.selected = {}
                self.snp_to_gene = {}
                for gene, rows in self.genes.iteritems():
                    for snp, row in rows.iteritems():
                        self.snp_to_gene[snp] = gene

            def __call__(self, i, comps):
                if i == 0:
                    return
                F = GTExSNPFile.SNPTF
                ref = comps[F.Ref_b37]
                eff = comps[F.Alt]
                snp = comps[F.RS_ID_dbSNP142_CHG37p13]
                if snp in self.snp_to_gene:
                    gene = self.snp_to_gene[snp]
                    row = self.genes[gene][snp]
                    if not gene in self.selected:
                        self.selected[gene] = {}
                    self.selected[gene][snp] = snp_fixed_row(row, snp, ref, eff)

        snp_callback = FixSNPCallback(selected)
        snp_iterator = Utilities.CSVFileIterator(gtex_snp_path, delimiter="\t", compressed=True)  # header in the file is just wrong., header=GTExSNPFile.SNPTF.HEADER, compressed=True)
        snp_iterator.iterate(snp_callback)
        selected = snp_callback.selected
        logging.info("%d survived after snp file", len(selected))
        del snp_callback

        logging.info("Fixing rows")
        fix_rows(selected, pvalues)


        logging.info("Saving entries")
        connection = Utilities.connect(db_output_path)
        Utilities.setup_db(connection)
        Utilities.insert_entries(connection, selected)