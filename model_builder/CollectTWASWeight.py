#! /usr/bin/env python

import os
import logging
import Logging
import Utilities
import GencodeFile
import TWASFormat

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

        sub_path = TWASFormat.build_subpaths(folder, content)

        map_path = sub_path + ".wgt.map"
        snps = TWASFormat.load_map(map_path)

        weights = TWASFormat.build_weights(sub_path)

        rows = []
        gene_id = gene_names[content]
        for i, snp in enumerate(snps):
            w = weights[i]
            row = (snp[TWASFormat.MTF.snp], gene_id, content, w, snp[TWASFormat.MTF.a1], snp[TWASFormat.MTF.a2] )
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