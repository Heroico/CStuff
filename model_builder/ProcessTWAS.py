#!/usr/bin/env python

import os
import logging
import shutil
import Utilities
import Logging
import TWASFormat

class GTF:
    snp=0
    a1=1
    a2=2
    z=3

class GiantBMICallback(object):
    class TF:
        SNP=0
        A1=1 #a1 is trait increasing
        A2=2
        Freq1_Hapmap=3
        b=4
        se=5
        p=6
        N=7

    def __init__(self):
        self.results = {}

    def __call__(self, i, comps):
        if i==0:
            return

        b = comps[GiantBMICallback.TF.b]
        se = comps[GiantBMICallback.TF.se]
        if b == "NA" or se == "NA" or float(se) == 0:
            return

        snp = comps[GiantBMICallback.TF.SNP]
        a1 = comps[GiantBMICallback.TF.A1]
        a2 = comps[GiantBMICallback.TF.A2]
        z = float(b)/float(se)

        self.results[snp] = (snp, a1, a2, z)

def build_twas_gene_folder(weight_folder, working_folder, gwas_results, content):
    print "yo"
    sub_path = TWASFormat.build_subpaths(weight_folder, content)

    map_path = sub_path + ".wgt.map"
    snps = TWASFormat.load_map(map_path)

    cor_path = sub_path + ".wgt.cor"
    cors = TWASFormat.load_cor(cor_path)

    ld_path = sub_path + ".wgt.ld"

    zscores = []
    for i,snp in enumerate(snps):
        zscore = 0
        rsid = snp[TWASFormat.MTF.snp]
        a1 = snp[TWASFormat.MTF.a1]
        a2 = snp[TWASFormat.MTF.a2]
        if rsid in gwas_results:
            result = gwas_results[rsid]
            if result[GTF.a1] == a1 and result[GTF.a2] == a2:
                zscore = result[GTF.z]
            elif result[GTF.a1] == a2 and result[GTF.a2] == a1:
                zscore =  -1.0*result[GTF.z]
            else:
                cors[i] = 0
        else:
            cors[i] = 0

        if zscore == 0:
            continue
        zscore = str(zscore)
        pos = snp[TWASFormat.MTF.pos]
        zscore_row = (rsid, pos, a1, a2, zscore)
        zscores.append(zscore_row)

    working_path = os.path.join(working_folder, content)
    os.makedirs(working_path)

    sub_working_path = os.path.join(working_path, content)

    working_cor = sub_working_path +".wgt.cor"
    with open(working_cor, "w") as file:
        for cor in cors:
            line = str(cor)+"\n"
            file.write(line)

    working_ld_path = sub_working_path + ".wgt.ld"
    shutil.copy(ld_path, working_ld_path)

    working_map_path = sub_working_path + ".wgt.map"
    shutil.copy(map_path, working_map_path)

    zscore_path = sub_working_path + ".zscore"
    with open(zscore_path, "w") as file:
        file.write("SNP_ID SNP_Pos Ref_Allele Alt_Allele Z-score\n")
        for z in zscores:
            line = " ".join(z) + "\n"
            file.write(line)


def parse_folder(input_folder, weight_folder, working_folder, gwas_results):
    contents = os.listdir(weight_folder)
    logging.info("processing folder")
    print len(contents)
    for content in contents:
        if content != "CCDC101":
            continue

        build_twas_gene_folder(weight_folder, working_folder, gwas_results, content)
        break


class RunTWAS(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        if os.path.exists(self.args.output):
            logging.info("%s already exists, delete it if you want ity done again", self.args.output)
            return

        gwas_parser = GiantBMICallback()
        logging.info("Processing gwas file")
        Utilities.parse_file(self.args.gwas_file, gwas_parser)
        gwas_results = gwas_parser.results
        parse_folder(self.args.input_folder, self.args.weight_folder, self.args.working_folder, gwas_results)

if __name__ == "__main__":
    class Args(object):
        def __init__(self):
            self.gencode_file = "data/gencode.v22.annotation.gtf.gz"
            self.input_folder = "/home/heroico/Documents/Projects/Chicago/3rd/TWAS"
            self.weight_folder = "/home/heroico/Documents/Projects/Chicago/3rd/TWAS/WEIGHTS_YFS"
            self.working_folder = "/home/heroico/Documents/Projects/Chicago/3rd/TWAS/WORKING"
            self.gwas_file = "/home/heroico/Documents/Projects/Chicago/MetaXcan/software/data/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
            self.output = "results/YFS_GIANT.csv"
            self.verbosity = 10

    args = Args()
    Logging.configureLogging(int(args.verbosity))

    work = RunTWAS(args)
    work.run()