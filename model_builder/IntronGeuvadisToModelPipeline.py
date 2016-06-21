#! /usr/bin/env python
import logging
import os

import Logging
import GEUVADISDosageFile
import GEUVADISExpressionFile
import LCEFileInfo

def run(args):
    logging.info("Loading geuvadis samples")
    geuvadis_sample_path = os.path.join(args.geuvadis_dosage_path,"samples.txt")
    geuvadis_samples = GEUVADISDosageFile.sample_ids(geuvadis_sample_path)

    logging.info("Loading intron samples")
    intron_samples = LCEFileInfo.sample_ids(args.intron_pheno_path)

    intersection = sorted([x for x in geuvadis_samples if x in intron_samples])
    logging.info("Found %d samples in common", len(intersection))

    logging.info("Building GEUVADIS dosage")
    GEUVADISDosageFile.dosage_to_pipeline_genotype(
            intersection, args.geuvadis_dosage_path, args.geuvadis_dosage_output)

    logging.info("Building GEUVADIS expression")
    GEUVADISExpressionFile.to_pipeline_expression(
            intersection, args.geuvadis_pheno_path, args.geuvadis_expression_output)

    logging.info("Building Leafcutter expression")
    LCEFileInfo.to_pipeline_expression(intersection,args.intron_pheno_path, args.intron_expression_output, args.intron_annotation_output)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Build Intron Quantification and Geuvadis phenotype to transcriptome model pipeline input")

    parser.add_argument("--geuvadis_dosage_path",
                        help="path to geuvadis dosage files",
                        default="data/dosagefiles-hapmap2")

    parser.add_argument("--geuvadis_pheno_path",
                        help="path to geuvadis file",
                        default="data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")

    parser.add_argument("--intron_pheno_path",
                        help="path to intron file",
                        default="data/GE_CEU.txt.sorted.gz")

    parser.add_argument("--gtex_snp_path",
                        help="path to GTEx snp file",
                        default="data/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz")

    parser.add_argument("--geuvadis_dosage_output",
                        help="output where geuvadis dosage will be saved",
                        default="results/geuvadis_dosage.txt.gz")

    parser.add_argument("--geuvadis_expression_output",
                        help="output where geuvadis expression will be saved",
                        default="results/geuvadis_expression.txt.gz")

    parser.add_argument("--intron_expression_output",
                        help="output where intron expression will be saved",
                        default="results/intron_expression.txt.gz")

    parser.add_argument("--intron_annotation_output",
                        help="output where intron annotation will be saved",
                        default="results/intron_annotation.txt.gz")

    parser.add_argument("--expect_throw",
                    help="Debug mode for incomplete gzipped file",
                    action="store_true",
                    default=False)


    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    run(args)