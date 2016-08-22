#! /usr/bin/env python

import GencodeFile
import Logging
import json
import logging
import os

def read_input_table(path):
    lines = []
    with open(path) as file:
        for l in file:
            lines.append(l.strip())
    return lines

def save_results(path, stats, input_lines):
    with open(path, "w") as file:
        keys_1 = stats.gene_type.keys()
        keys_2 = stats.transcript_type.keys()
        print keys_1
        print keys_2

        for l in input_lines:
            comps = l.split("|")
            k = comps[0].strip() if len(comps[0]) else comps[1].strip()
            # print comps
            print "-"+k+"-"
            g = stats.gene_type[k]["count"] if k in stats.gene_type else "NA"
            t = stats.transcript_type[k]["count"] if k in stats.transcript_type else "NA"
            line = l + " | " + str(g) +" | "+ str(t) + "|\n"
            file.write(line)

        #json.dump({"gene_types":callback.gene_types, "transcript_types":callback.transcript_type}, file, sort_keys=True, indent=2)
def run(args):
    class Callback(object):
        def __init__(self):
            self.gene_type = {}
            self.transcript_type = {}

        def __call__(self, gencode):
            if not gencode.gene_type in self.gene_type:
                self.gene_type[gencode.gene_type] = {"count":0}
            self.gene_type[gencode.gene_type]["count"] = self.gene_type[gencode.gene_type]["count"] + 1

            if not gencode.transcript_type in self.transcript_type:
                self.transcript_type[gencode.transcript_type] = {"count":0}
            self.transcript_type[gencode.transcript_type]["count"] = self.transcript_type[gencode.transcript_type]["count"] + 1

    if os.path.exists(args.results_file):
        logging.info("Output file already exists.")
        exit(0)

    input_lines = read_input_table(args.input_table)

    callback = Callback()
    GencodeFile.parse_gencode_file(args.gencode_file, callback, only_genes=False)

    folder = os.path.split(args.results_file)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)

    save_results(args.results_file, callback, input_lines)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        "Build Intron Quantification and Geuvadis phenotype to transcriptome model pipeline input")


    parser.add_argument("--gencode_file",
                    help="output where intron annotation will be saved",
                    default="data/gencode.v19.annotation.gtf.gz")

    parser.add_argument("--input_table",
                    help="table to be filled with stats",
                    default="data/input.t")

    parser.add_argument("--results_file",
                    help="output where intron annotation will be saved",
                    default="results/gencode.stats.t")

    parser.add_argument("--verbosity",
                    help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                    default="10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    logging.info("running")
    run(args)
    logging.info("ran")