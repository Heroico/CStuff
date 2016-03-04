#! /usr/bin/env python

import os
import logging
import math
import scipy.stats as stats
import PB8KFileInfo
import VarianceFile
import Logging
import Utilities

def row_from_comps(comps, TF, var, N):
    snp = comps[TF.SNPName]
    if not snp in var:
        return None
    chr = "chr"+comps[TF.SNPChr]
    pos = comps[TF.SNPChrPos]
    reference_allele, effect_allele = PB8KFileInfo.alleles(TF, comps)
    zscore = float(comps[TF.OverallZScore])
    se = math.sqrt(float(var[snp])/N)
    beta = zscore*se
    OR = str(math.exp(beta))
    #check ldpred coord_genotype parse_sum_stats_standard() for what is a1 and a2
    p = str(stats.norm.sf(abs(zscore)) * 2)
    row = (chr, snp, reference_allele, effect_allele, pos, OR, p,)
    return row


#hg19chrc    snpid    a1    a2    bp    or    p
def gather_pb8k(file_path, variances, N):
    class ParseCallback(object):
        def __init__(self, variance, N):
            self.variance = variance
            self.N = N
            self.index = {}
            self.results = {}

        def __call__(self, line_number, comps):
            if line_number==0:
                return
            gene = comps[PB8KFileInfo.TF1.HUGO]
            if "," in gene:
                multiple_genes = gene.split(",")
                for g in multiple_genes:
                    row = row_from_comps(comps, PB8KFileInfo.TF1, self.variance, self.N)
                    self.process_row(g, row)
            else:
                row = row_from_comps(comps, PB8KFileInfo.TF1, self.variance, self.N)
                self.process_row(gene, row)

        def process_row(self, gene, row):
            if row is None:
                return
            snp = row[1]
            self.index[snp] = gene

            if not gene in self.results:
                self.results[gene] = []
            self.results[gene].append(row)

    callback = ParseCallback(variances, N)
    Utilities.parse_file(file_path, callback, expect_throw=True)
    return callback.index, callback.results

def save_slice(output_path, gene, rows):
    p = os.path.join(output_path, gene)
    os.makedirs(p)
    p = os.path.join(p, gene+".txt")
    with open(p, "w") as file:
        file.write("hg19chrc snpid a1 a2 bp or p\n")
        for r in rows:
            line = " ".join(r)+"\n"
            file.write(line)

class Slice(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        if os.path.exists(self.args.output_folder):
            logging.info("Output folder already exists, delete it if you want it done again")
            return

        os.makedirs(self.args.output_folder)

        logging.info("Loading variance")
        variances = VarianceFile.load_variance(self.args.variance_path)

        logging.info("Processing GWAS")
        index, results = gather_pb8k(self.args.pb8k_path, variances, self.args.sample_size)

        logging.info("Writing slices")
        for gene, rows in results.iteritems():
            save_slice(self.args.output_folder, gene, rows)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Slice PB8K GWAS and Thousand Genomes files for ldpred input")

    parser.add_argument("--pb8k_path",
                        help="path to PB8K file",
                        default="data/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt.gz")

    parser.add_argument("--variance_path",
                        help="path to Variance file",
                        default="data/VAR_TGF_EUR_PB8K.txt.gz")

    parser.add_argument("--sample_size",
                        help="Number of individuals",
                        type=int,
                        default=5311)

    parser.add_argument("--thousand_genomes_folder",
                        help="path to thousand genomes phase 3 impute2 format files",
                        default="/home/heroico/Documents/Projects/Chicago/Data/1000GP_Phase3")

    parser.add_argument("--output_folder",
                        help="path where results will be saved",
                        default="results/sliced")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = Slice(args)
    work.run()