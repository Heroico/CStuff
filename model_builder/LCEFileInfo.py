import gzip
import logging
import os
import cStringIO

#Leafcutter expression table format
class LEFTF(object):
    CHR=0
    START=1
    END=2
    ID=3
    FIRST_INDIVIDUAL=4

#id to sample index, column in this case
def sample_ids(intron_path):
    samples = []
    with gzip.open(intron_path) as file:
        first = file.readline().strip().split()
        samples = first[LEFTF.FIRST_INDIVIDUAL:]
    return samples

class PipelineExpressionOutput(object):
    def __init__(self, output_file, samples, selected):
        self.output_file = output_file
        self.samples = samples
        self.selected = selected
        self.selected_index = {x:i for i,x in enumerate(samples) if x in selected}
        logging.info("%s samples found at Leafcutter expression pipeline output", len(self.selected))

    def __call__(self, line):
        id = line[LEFTF.ID]
        comps = id.split(":")
        intron_id = "_".join(comps[1:])
        dosage = line[LEFTF.FIRST_INDIVIDUAL]
        dosage = [dosage[self.selected_index[x]] for x in self.selected]
        output = "\t".join([intron_id]+dosage)+"\n"
        self.output_file.write(output)

class PipelineGeneAnnotationOutput(object):
    def __init__(self, output_file):
        self.output_file = output_file

    def __call__(self, line):
        chr = "chr"+line[LEFTF.CHR]
        start_pos = line[LEFTF.START]
        end_pos = line[LEFTF.END]
        id = line[LEFTF.ID]
        comps = id.split(":")
        name = "_".join(comps[1:])

#chr1    NA  gene      11869   14409   .       +       .       gene_id "668593.672093.clu_46640";
#  transcript_id "668593.672093.clu_46640"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "668593.672093.clu_46640";
#  transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "668593.672093.clu_46640";
#  level 2; tag "basic"; transcript_support_level "1"; havana_gene "668593.672093.clu_46640"; havana_transcript "668593.672093.clu_46640";
        l = cStringIO.StringIO()
        l.write(chr), l.write("\t")
        l.write("NA\tgene\t")
        l.write(start_pos), l.write("\t")
        l.write(end_pos), l.write("\t")
        l.write(".\t+\t.\t")
        l.write('gene_id "'), l.write(name), l.write('";\t')
        l.write('transcript_id "'), l.write(name), l.write('";\t')
        l.write('gene_type "protein_coding";\t')
        l.write('gene_status "KNOWN";\t')
        l.write('gene_name "'), l.write(name), l.write('";\t')
        l.write('transcript_type "processed_transcript";\t')
        l.write('gene_status "KNOWN";\t')
        l.write('transcript_name "'), l.write(name), l.write('";\t')
        l.write('level 2;\ttag "basic";\ttranscript_support_level "1";\t')
        l.write('havana_gene "'), l.write(name), l.write('";\t')
        l.write('havana_transcript "'), l.write(name), l.write('";\n')
        self.output_file.write(l.getvalue())

class LeafcutterExpressionIterator(object):
    def __init__(self, path):
        self.path = path

    def iterate(self, callbacks):
        with gzip.open(self.path) as file:
            for i,line in enumerate(file):
                if i==0:
                    continue
                comps = line.strip().split()
                chr = comps[LEFTF.CHR]
                start = comps[LEFTF.START]
                end = comps[LEFTF.END]
                id = comps[LEFTF.ID]
                dosage = comps[LEFTF.FIRST_INDIVIDUAL:]
                data = (chr,start,end,id,dosage)

                for c in callbacks:
                    c(data)

def to_pipeline_expression(selected, file, expression_output_path, annotation_ouput_path):
    do_expression = True
    if os.path.exists(expression_output_path):
        logging.info("%s: intron expression output exists, delete it if you want it done again", expression_output_path)
        do_expression = False

    do_annotation = True
    if os.path.exists(annotation_ouput_path):
        logging.info("%s: intron annotation output exists, delete it if you want it done again", annotation_ouput_path)
        do_annotation = False

    if not do_expression and not do_annotation:
        return

    folder = os.path.split(expression_output_path)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    samples = sample_ids(file)
    callbacks = []
    if do_expression:
        expression_file = gzip.open(expression_output_path, "wb")
        e_header = ["Id"]+selected
        e_header = "\t".join(e_header)+"\n"
        expression_file.write(e_header)
        expression_callback = PipelineExpressionOutput(expression_file, samples, selected)
        callbacks.append(expression_callback)

    if do_annotation:
        annotation_file = gzip.open(annotation_ouput_path, "wb")
        annotation_file.write("##description: evidence-based annotation of the human genome (GRCh38), version 22 (Ensembl 79)\n")
        annotation_file.write("##provider: GENCODE\n")
        annotation_file.write("##contact: gencode@sanger.ac.uk\n")
        annotation_file.write("##format: gtf\n")
        annotation_file.write("##date: 2015-03-06\n")
        annotation_callback = PipelineGeneAnnotationOutput(annotation_file)
        callbacks.append(annotation_callback)

    iterator = LeafcutterExpressionIterator(file)
    iterator.iterate(callbacks)

    if do_expression:
        expression_file.close()

    if do_annotation:
        annotation_file.close()
