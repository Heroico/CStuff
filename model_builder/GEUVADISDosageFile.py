import gzip
import os
import logging
import Utilities

#id to sample index, row in this case
def sample_ids(geuvadis_sample_path):
    samples = []
    with open(geuvadis_sample_path) as file:
        for i,l in enumerate(file):
            comps = l.strip().split()
            samples.append(comps[0])
    return samples

# Geuvadis dosgae file table format
class GDFTF(object):
    CHR=0
    ID=1
    POS=2
    REFERENCE_ALLELE=3
    EFFECT_ALLELE=4
    AVERAGE=5
    FIRST_DOSAGE=6

class PipelineOutputCollector(object):
    def __init__(self, output_file, samples, selected):
        self.output_file = output_file
        self.samples = samples
        self.selected = selected
        self.selected_index = {x:i for i,x in enumerate(samples) if x in selected}
        logging.info("%s samples found at GEUVADIS dosage pipeline output", len(self.selected))

    def __call__(self, line):
        #chr1 rs10399749 55299 C T 0.164302600472813 ...
        chr = line[GDFTF.CHR]
        chr_number = chr.split("chr")[1]
        pos = line[GDFTF.POS]
        ref = line[GDFTF.REFERENCE_ALLELE]
        eff = line[GDFTF.EFFECT_ALLELE]
        dosage = line[GDFTF.FIRST_DOSAGE]
        dosage = [dosage[self.selected_index[x]] for x in self.selected]
        # chr_pos_ref_eff_build
        id = "_".join([chr_number, pos, ref, eff, "b37"])
        output = "\t".join([id]+dosage)+"\n"
        self.output_file.write(output)

class GEUVADISDosageFileIterator(object):
    def __init__(self, path):
        self.path = path

    def iterate(self, callback):
        with gzip.open(self.path) as dosage_file:
            for i,l in enumerate(dosage_file):
                comps = l.strip().split()
                #chr1 rs10399749 55299 C T 0.164302600472813 ...
                chr = comps[GDFTF.CHR]
                id = comps[GDFTF.ID]
                pos = comps[GDFTF.POS]
                ref = comps[GDFTF.REFERENCE_ALLELE]
                eff = comps[GDFTF.EFFECT_ALLELE]
                average = comps[GDFTF.AVERAGE]
                dosage = comps[GDFTF.FIRST_DOSAGE:]
                line = (chr, id, pos, ref, eff, average, dosage)
                callback(line)

def dosage_to_pipeline_genotype(selected, geuvadis_folder_path, output_path):
    if os.path.exists(output_path):
        logging.info("GEUVADIS dosage output exists, delete it if you want it done again")
        return

    folder = os.path.split(output_path)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    geuvadis_sample_path = os.path.join(geuvadis_folder_path,"samples.txt")
    file_samples = sample_ids(geuvadis_sample_path)
    with gzip.open(output_path, "wb") as output_file:
        header = ["Id"]+[x for i,x in enumerate(selected)]
        header = "\t".join(header)+"\n"
        output_file.write(header)
        contents = os.listdir(geuvadis_folder_path)
        for c in contents:
            if "samples" in c:
                continue
            path = os.path.join(geuvadis_folder_path, c)
            logging.info("Processing %s", path)
            callback = PipelineOutputCollector(output_file, file_samples, selected)
            iterator = GEUVADISDosageFileIterator(path)
            iterator.iterate(callback)
