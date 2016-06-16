import logging
import gzip
import os

#Geuvadis expression file table format
class GEFTF(object):
    TARGET_ID=0
    GENE_SYMBOL=1
    CHR=2
    COORD=3
    FIRST_INDIVIDUAL=4

class PipelineOutput(object):
    def __init__(self, output_file, samples, selected):
        self.output_file = output_file
        self.samples = samples
        self.selected = selected
        self.selected_index = {x:i for i,x in enumerate(samples) if x in selected}
        logging.info("%s samples found at GEUVADIS expression pipeline output", len(self.selected))
    
    def __call__(self, line):
        gene_id = line[GEFTF.TARGET_ID]
        dosage = line[GEFTF.FIRST_INDIVIDUAL]
        dosage = [dosage[self.selected_index[x]] for x in self.selected]
        output = "\t".join([gene_id]+dosage)+"\n"
        self.output_file.write(output)

class GEUVADISExpressionIterator(object):
    def __init__(self, path):
        self.path = path
        self.sample_ids = None

    def iterate(self, calback):
        with gzip.open(self.path) as file:
            for i,line in enumerate(file):
                if i==0:
                    continue
                comps = line.strip().split()
                target_id = comps[GEFTF.TARGET_ID]
                gene_symbol = comps[GEFTF.GENE_SYMBOL]
                chr = comps[GEFTF.CHR]
                coord = comps[GEFTF.COORD]
                dosage = comps[GEFTF.FIRST_INDIVIDUAL:]
                data = (target_id,gene_symbol,chr,coord,dosage)
                calback(data)

def sample_ids(path):
    samples = []
    with gzip.open(path) as file:
        comps = file.readline().strip().split()
        samples = comps[GEFTF.FIRST_INDIVIDUAL:]
    return samples

def to_pipeline_expression(selected, file, output_path):
    if os.path.exists(output_path):
        logging.info("GEUVADIS expression output exists, delete it if you want it done again")
        return

    folder = os.path.split(output_path)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    samples = sample_ids(file)
    with gzip.open(output_path, "wb") as output_file:
        header = ["Id"]+selected
        header = "\t".join(header)+"\n"
        output_file.write(header)
        callback = PipelineOutput(output_file, samples, selected)
        iterator = GEUVADISExpressionIterator(file)
        iterator.iterate(callback)

    

