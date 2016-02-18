#! /usr/bin/env python
import sqlite3
import logging
import os
import gzip
import math
import Logging
import GencodeFile
import Utilities

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

def row_from_comps(gene, comps, TF):
    snp = comps[TF.SNPName]
    alleles = {a for a in comps[TF.SNPType].split("/")}

    #TODO: figure this out
    assessed = comps[TF.AlleleAssessed]
    for a in alleles:
        if a == assessed:
            effect_allele = a
        else:
            reference_allele = a
    zscore = comps[TF.OverallZScore]
    row = (snp, gene, zscore, reference_allele, effect_allele,)
    return row

#add ensemble id. "gene" becomes a "gene name", conceptually
def process_row(gene, row, genes, only_best_snp =None):
    if not gene in genes:
        genes[gene] = []
    rows = genes[gene]
    if only_best_snp:
        if len(rows):
            r = rows[0]
            if math.fabs(float(r[2])) < math.fabs(float(row[2])):
                del rows[0]
                rows.append(row)
        else:
            rows.append(row)
    else:
        rows.append(row)

def parse_input_file(TF, connection, input_file, gencode_file, fdr_filter=None, use_variance=None, sample_size=None, only_best_snp=None):
    genes = {}
    logging.info("Opening pheno phile")
    with gzip.open(input_file) as file:
        for i,line in enumerate(file):
            if i==0:
                continue

            comps = line.strip().split()

            if fdr_filter:
                fdr = comps[TF.FDR]
                if float(fdr) > fdr_filter:
                    snp = comps[TF.SNPName]
                    logging.log(9,"Snp %s doesn't pass fdr filter: %s", snp, fdr)
                    continue

            gene = comps[TF.HUGO]
            if "," in gene:
                multiple_genes = gene.split(",")
                for g in multiple_genes:
                    row = row_from_comps(g, comps, TF)
                    process_row(g, row, genes, only_best_snp)
            else:
                row = row_from_comps(gene, comps, TF)
                process_row(gene, row, genes, only_best_snp)

    logging.info("Opening gencode file")

    class GenCodeCallback(object):
        def __init__(self, genes):
            self.genes = genes
            self.selected = {}

        def __call__(self, gencode):
            if gencode.name in self.genes:
                rows = self.genes[gencode.name]
                self.selected[gencode.name] = [(row[0], gencode.ensemble_version, row[1], row[2], row[3], row[4]) for row in rows]
    callback = GenCodeCallback(genes)
    GencodeFile.parse_gencode_file(gencode_file, callback)
    genes = callback.selected

    if use_variance:
        logging.info("Opening variance file")
        vars = {}
        with gzip.open(use_variance, "rb") as var_file:
            for line in var_file:
                comps = line.strip().split(",")
                vars[comps[0]] = float(comps[1])
        keys = genes.keys()
        for key in keys:
            rows = genes[key]
            new_rows = []
            for r in rows:
                snp = r[0]
                if not snp in vars:
                    continue
                v = vars[snp]
                std = math.sqrt(v/sample_size)
                new_rows.append([r[0], r[1], r[2], str(float(r[3])*std), r[4], r[5]])
            genes[key] = new_rows

    Utilities.insert_entries(connection, genes)

class BuildModel(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        if os.path.exists(self.args.output_file):
            logging.info("DB already there, delete it if you want it done again")
            return

        connection = Utilities.connect(self.args.output_file)

        Utilities.setup_db(connection)
        parse_input_file(TF1, connection, self.args.input_file, self.args.gencode_file, self.args.fdr_filter,
                         self.args.use_variance, self.args.sample_size, self.args.only_best_snp)

if __name__ == "__main__":
    class Args(object):
        def __init__(self):
            self.input_file = "data/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt.gz"
            self.gencode_file = "data/gencode.v22.annotation.gtf.gz"
            self.output_file = "results/PB8K_beta_f_best_snp.db"
            self.fdr_filter = float(0.05)
            self.use_variance = "data/VAR_TGF_EUR_PB8K.txt.gz"
            self.sample_size = 5311
            self.only_best_snp = True
            self.verbosity = 10

    args = Args()
    Logging.configureLogging(int(args.verbosity))

    work = BuildModel(args)
    work.run()