#! /usr/bin/env python
__author__ = 'heroico'

import os
import logging
import time
import sys
from subprocess import call


def configureLogging(level=5, target=sys.stderr):
    logger = logging.getLogger()
    logger.setLevel(level)

    # create console handler and set level to info
    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def run_job(input_path, output_path, job_script_path, logs_folder):
    if os.path.exists(job_script_path):
        logging.info("File %s already exisits, skipping", job_script_path)

    with open(job_script_path, "w") as file:
        job = "#!/bin/bash" + "\n\n"
        job += "#PBS -N " + job_script_path + "\n\n"
        job += "#PBS -S /bin/bash" + "\n\n"
        job += "#PBS -l walltime=1:00:00" + "\n\n"
        job += "#PBS -l nodes=1:ppn=1" + "\n\n"
        job += "#PBS -l mem=4gb" + "\n\n"
        job += "#PBS -o " + logs_folder + "/${PBS_JOBNAME}.o${PBS_JOBID}.log" + "\n\n"
        job += "#PBS -e " + logs_folder + "/${PBS_JOBNAME}.e${PBS_JOBID}.err" + "\n\n"
        job += "cd $PBS_O_WORKDIR" + "\n\n"
        job += "module load plink/1.90" +"\n\n"
        job += "plink --vcf " + input_path + " --out " + output_path
        file.write(job)
    retry = 0

    # while call(["qsub",job_script_path]) != 0 and retry < 10:
    #     logging.info("retry %s %i-th attempt", job_script_path, retry)
    #     retry += 1
    #     time.sleep(0.1)
    # time.sleep(0.1)

def run(args):
    if not os.path.exists(args.intermediate_folder):
        os.makedirs(args.intermediate_folder)
        logs_folder = os.path.join(args.intermediate_folder, "logs")
        os.makedirs(logs_folder)
    else:
        logging.info("Intermediate folder exists, delete it you want this run again")
        return

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    else:
        logging.info("Output folder exists, delete it you want this run again")
        return

    contents = os.listdir(args.input_folder)
    contents = [x for x in contents if "vcf.gz" in x and not ".tbi" in x]
    i = 0
    for c in contents:
        input_path = os.path.join(args.input_folder, c)

        output_name = c.split(".vcf.gz")[0]
        output_path = os.path.join(args.output_folder, output_name)

        job_name = "job_%i.sh" % (i,)
        job_path = os.path.join(args.intermediate_folder, job_name)

        run_job(input_path, output_path, job_path, logs_folder)
        i = i+1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Run plink on 1000G files")

    parser.add_argument("--input_folder",
                        help="folder with 1000G data",
                        default="/group/im-lab/nas40t2/abarbeira/data/1000G")

    parser.add_argument("--output_folder",
                        help="folder with converted results",
                        default="/group/im-lab/nas40t2/abarbeira/data/P_1000G/results")

    parser.add_argument("--intermediate_folder",
                        help="folder with scratch",
                        default="/group/im-lab/nas40t2/abarbeira/data/P_1000G/intermediate")


    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    configureLogging(int(args.verbosity))

    run(args)