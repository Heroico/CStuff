#!/usr/bin/env python
import os
import sys
import logging
import time
from subprocess import call, Popen, PIPE

def configureLogging(level=5, log_file=None, target=sys.stderr):
    logger = logging.getLogger()
    logger.setLevel(level)

    # create console handler and set level to info
    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if log_file:
        if os.path.exists(log_file):
            os.remove(log_file)
        fh = logging.FileHandler(log_file)
        fh.setLevel(level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

def build_command(args, weight):
    name = os.path.split(weight)[1]
    name = name.split(".db")[0]
    output = os.path.join(args.output_folder, name) + ".expr.txt"

    command = "{} \\\n".format(args.predict_expression_command)
    command += "--dosages {} \\\n".format(args.dosages_folder)
    command += "--weights {} \\\n".format(weight)
    if args.gtex_snps:
        command += "--gtex_snps {} \\\n".format(args.gtex_snps)
        if not args.dosages_prefix:
            prefix = [x.split(".gz")[0] for x in os.listdir(args.dosages_folder)
                      if ("_Analysis.snps.txt.gz" in x and x.split("_Analysis.snps.txt.gz")[0] in name)][0]
            command += "--dosages_prefix {} \\\n".format(prefix)

    if args.dosages_prefix:
        command += "--dosages_prefix {} \\\n".format(args.dosages_prefix)

    command += "--output {}".format(output)

    return name, command

def job_header(job_name, logs_folder="logs", partition="broadwl", time="2:00:00", mem="4096", ntasks="1", email=None):

    header = "#!/bin/bash \n"
    header += "#SBATCH --job-name={} \n".format(job_name)

    if logs_folder:
        log_name = os.path.join(logs_folder, job_name)+".log"
        err_name = os.path.join(logs_folder, job_name) + ".err"

        header += "#SBATCH --output={} \n".format(log_name)
        header += "#SBATCH --error={} \n".format(err_name)

    if email:
        header += "#SBATCH --mail-user={} \n".format(email)
        header += "#SBATCH --mail-type=END,FAIL \n".format(email)
# SBATCH --array=1-22
    header+="#SBATCH --partition={} \n".format(partition)
    header+="#SBATCH --mem-per-cpu={} \n".format(mem)
    header+="#SBATCH --ntasks={} \n".format(ntasks)
    header+="#SBATCH --time={} \n".format(time)
    return header

def submit_command(job_folder, log_folder, job_name, command, fake=False, serialize_local=False, partition=None):
    path = os.path.join(job_folder, job_name+ ".sbatch")

    if os.path.exists(path):
        logging.info("job %s already exists, we assume it to be queued in the past", job_name)
        return None, None

    c = job_header(job_name,log_folder, time="12:00:00", mem=4096, partition=partition, email="mundoconspam@gmail.com")
    c += "\n"
    c += "module load python\n\n"
    c += command

    with open(path, "w") as file:
        file.write(c)

    if fake: return job_name, None

    if serialize_local:
        logging.info("Running %s", path)
        call(["bash", path])
        return job_name

    return path, submit_job(path, job_name=job_name)

def submit_job(path ,job_name=None):
    job_name = os.path.split(path)[1] if not job_name else job_name
    retry = 0
    command = ["sbatch",path]
    def submit():
        logging.log(8,"Submitting Command: %s", " ".join(command))
        proc = Popen(command,stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
        return exitcode, out, err

    exitcode, out, err = submit()
    while exitcode != 0 and retry < 10:
        logging.info("retry %s %i-th attempt", job_name, retry)
        exitcode, out, err = submit()
        retry += 1
        time.sleep(0.5)
    time.sleep(0.5)
    r=None
    if exitcode==0:
        r=out.strip()
        r=str(r)
        r=r.split("Submitted batch job ")[1]

    return r

def _args():
    class Args1(object):
        def __init__(self):
            self.logs = "logs"
            self.jobs = "jobs"
            self.weights_folder = "/project/haky/im-lab/nas40t2/abarbeira/predictdb/GTEx-V6p-HapMap-2016-09-08"
            self.predict_expression_command=  "/project/haky/im-lab/nas40t2/abarbeira/software/PrediXcan/Software/predict_gene_expression.py"
            self.dosages_folder = "/project/haky/im-lab/nas40t2/kaanan/BD"
            self.dosages_prefix = "BD_chr"
            self.gtex_snps = None
            self.partition="broadwl"
            self.output_folder = "results_bd"
            # self.dosages_folder = "/project/haky/im-lab/nas40t2/kaanan/T1D"
            # self.dosages_prefix = "T1D_chr"
            # self.partition="broadwl"
            # self.output_folder = "results_t1d"
            self.verbosity = 10
            self.fake_submission = False

    class Args2(object):
        def __init__(self):
            self.logs = "logs"
            self.jobs = "jobs"
            self.weights_folder = "/project/haky/im-lab/nas40t2/abarbeira/predictdb/GTEx-V6p-HapMap-2016-09-08"
            self.predict_expression_command=  "/project/haky/im-lab/nas40t2/abarbeira/software/PrediXcan/Software/predict_gene_expression.py"
            self.dosages_folder = "/project/haky/Data/GTEx/V6p/genotypes/"
            self.dosages_prefix = "Muscle_Skeletal"
            self.gtex_snps = "/project/haky/Data/GTEx/V6p/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt.gz"
            self.partition="broadwl"
            self.output_folder = "results_gtex"
            self.verbosity = 10
            self.fake_submission = False

    return Args2()

def run(args):
    if not os.path.exists(args.logs): os.makedirs(args.logs)
    if not os.path.exists(args.jobs): os.makedirs(args.jobs)
    if not os.path.exists(args.output_folder): os.makedirs(args.output_folder)

    weights = sorted([os.path.join(args.weights_folder,x) for x in os.listdir(args.weights_folder) if ".db" in x])
    for weight in weights:
        name, command = build_command(args, weight)
        path, status = submit_command(args.jobs, args.logs, name, command, fake = args.fake_submission, partition=args.partition)
        logging.info("%s|%s", path, str(status))

if __name__ == "__main__":
    args = _args()
    configureLogging(args.verbosity)
    run(args)
