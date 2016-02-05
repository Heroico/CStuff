#! /usr/bin/env python
import sqlite3
import logging
import os
import gzip
import Logging
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

def build_db(db_path):
    return sqlite3.connect(db_path)

def setup_db(connection):
    cursor = connection.cursor()
    cursor.execute("CREATE TABLE extra (gene TEXT, genename TEXT, R2 DOUBLE,  `n.snps` INTEGER)")
    cursor.execute("CREATE INDEX extra_gene ON extra (gene)")
    cursor.execute("CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
    cursor.execute("CREATE INDEX weights_rsid ON weights (rsid)")
    cursor.execute("CREATE INDEX weights_gene ON weights (gene)")
    cursor.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
    connection.commit()

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
def process_row(gene, row, genes):
    if not gene in genes:
        genes[gene] = []
    genes[gene].append(row)

def parse_input_file(connection, input_file, gencode_file, TF):
    genes = {}
    logging.info("Opening pheno phile")
    with gzip.open(input_file) as file:
        for i,line in enumerate(file):
            if i==0:
                continue

            comps = line.strip().split()

            gene = comps[TF.HUGO]
            if "," in gene:
                multiple_genes = gene.split(",")
                for g in multiple_genes:
                    row = row_from_comps(g, comps, TF)
                    process_row(g, row, genes)
            else:
                row = row_from_comps(gene, comps, TF)
                process_row(gene, row, genes)

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

    cursor = connection.cursor()

    logging.info("Inserting snp entries")
    data = []
    for rows in genes.values():
        cursor.executemany("INSERT INTO weights VALUES(?, ?, ?, ?, ?, NULL, NULL, NULL)", [(r[0], r[2], r[3], r[4], r[5]) for r in rows])

    logging.info("Inserting gene entries")
    for gene, rows in genes.iteritems():
        r = rows[0]
        cursor.execute("INSERT INTO extra VALUES(?, ?, ?, ?)", (r[1], r[2], "NA", len(r)))
    connection.commit()

class BuildModel(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        if os.path.exists(self.args.output_file):
            logging.info("DB already there, delete it if you want it done again")
            return

        connection = build_db(self.args.output_file)

        setup_db(connection)
        parse_input_file(connection, self.args.input_file, self.args.gencode_file, TF1)


if __name__ == "__main__":
    class Args(object):
        def __init__(self):
            self.input_file = "data/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt.gz"
            self.gencode_file = "data/gencode.v22.annotation.gtf.gz"
            self.output_file = "results/PB8K.db"
            self.verbosity = 10

    args = Args()
    Logging.configureLogging(int(args.verbosity))

    work = BuildModel(args)
    work.run()