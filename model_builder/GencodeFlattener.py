#! /usr/bin/env python

import GencodeFile
import Logging
import gzip
import logging
import os


def run(args):
    class Callback(object):
        def __init__(self, file):
            self.file = file
            header = "\t".join(["chr", "gene","gene_name","start_location","end_location"]) + "\n"
            file.write(header)

        def __call__(self, gencode):
            comps = [gencode.chromosome_name, gencode.gene_id, gencode.gene_name, gencode.start_location, gencode.end_location]
            line = "\t".join(comps) + "\n"
            file.write(line)

    if os.path.exists(args.results_file):
        logging.info("Output file already exists.")
        exit(0)

    folder = os.path.split(args.results_file)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)

    with gzip.open(args.results_file, "wb") as file:
        callback = Callback(file)
        GencodeFile.parse_gencode_file(args.gencode_file, callback, only_genes=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        "Build Intron Quantification and Geuvadis phenotype to transcriptome model pipeline input")


    parser.add_argument("--gencode_file",
                    help="output where intron annotation will be saved",
                    default="data/gencode.v19.annotation.gtf.gz")

    parser.add_argument("--results_file",
                    help="output where intron annotation will be saved",
                    default="results/gencode.v19.annotation.gtf.txt")

    parser.add_argument("--verbosity",
                    help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                    default="10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    logging.info("running")
    run(args)
    logging.info("ran")