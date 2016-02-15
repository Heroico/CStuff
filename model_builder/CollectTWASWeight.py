#! /usr/bin/env python

import os
import logging
import numpy
import Logging
import Utilities
import GencodeFile

def load_map(path):
    snps = []
    with open(path) as file:
        for i, line in enumerate(file):
            if i==0:
                continue
            comps = line.strip().split()
            row = (comps[0], comps[2], comps[3], )
            snps.append(row)
    return snps

def load_cor(path):
    cors = []
    with open(path) as file:
        for line in file:
            c = line.strip()
            cors.append(float(c))
    return cors

def load_ld(path):
    rows = []
    with open(path) as file:
        for line in file:
            row = [float(x) for x in line.strip().split()]
            rows.append(row)
    array = numpy.array(rows)
    i, j = numpy.indices(array.shape)
    array[i == j] += 0.01
    return array

def build_weights(sub_path):
    cor_path = sub_path + ".wgt.cor"
    cors = load_cor(cor_path)

    ld_path = sub_path + ".wgt.ld"
    ld = load_ld(ld_path)

    weights = calculate_weights(cors, ld)
    return weights

def calculate_weights(cors, ld):
    inv = numpy.linalg.inv(ld)
    dot = numpy.dot(cors, inv)
    return dot

def parse_folder(folder, db, gencode_file):

    contents = os.listdir(folder)

    logging.info("Processing gencode file")
    class GencodeCallback(object):
        def __init__(self, contents):
            self.contents = {gene:True for gene in contents}
            self.genes = {}
        def __call__(self, gencode):
            if gencode.name in self.contents:
                self.genes[gencode.name] = gencode.ensemble_version
    callback = GencodeCallback(contents)
    GencodeFile.parse_gencode_file(gencode_file, callback)
    gene_names = callback.genes

    logging.info("processing folder")
    genes = {}
    for content in contents:
        if not content in gene_names:
            logging.log(9, "Gene %s not in gencode", content)
            continue

        sub_path = os.path.join(folder, content)
        sub_path = os.path.join(sub_path, content)

        map_path = sub_path + ".wgt.map"
        snps = load_map(map_path)

        weights = build_weights(sub_path)

        rows = []
        gene_id = gene_names[content]
        for i, snp in enumerate(snps):
            w = weights[i]
            row = (snp[0], gene_id, content, snp[1], snp[2], w, )
            rows.append(row)
        genes[content] = rows

    Utilities.insert_entries(db, genes)

class Collect(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        if os.path.exists(self.args.output):
            logging.info("%s already exists, delete it if you want ity done again", self.args.output)
            return

        db = Utilities.connect(self.args.output)

        Utilities.setup_db(db)

        parse_folder(self.args.input_folder, db, self.args.gencode_file)

if __name__ == "__main__":
    class Args(object):
        def __init__(self):
            self.gencode_file = "data/gencode.v22.annotation.gtf.gz"
            self.input_folder = "/home/heroico/Documents/Projects/Chicago/3rd/TWAS/WEIGHTS_YFS"
            self.output = "results/YFS.db"
            self.verbosity = 10

    args = Args()
    Logging.configureLogging(int(args.verbosity))

    work = Collect(args)
    work.run()