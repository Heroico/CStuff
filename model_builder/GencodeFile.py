__author__ = 'heroico'
#trimmed from PredictDBAnalysis/gencode_input

import Utilities

K_NOT_GENES = ["transcript","exon","CDS","UTR","start_codon","stop_codon","Selenocysteine"];

# look at gencode http://www.gencodegenes.org/data_format.html
class GFTF:
    """gencode file table format"""
    CHROMOSOME = 0
    ANNOTATION_SOURCE = 1
    FEATURE_TYPE = 2
    N_START_LOCATION = 3
    N_END_LOCATION = 4
    SCORE = 5
    GENETIC_STRAND = 6
    GENOMIC_PHASE = 7
    KEY_VALUE_PAIRS = 8

    #there are several other key-value pairs but we are concerned with these
    GENE_ID = "gene_id"
    GENE_NAME = "gene_name"

class GFTFS:
    """gencode short file table format"""
    CHROMOSOME = 0
    FEATURE_TYPE = 1
    N_START_LOCATION = 2
    N_END_LOCATION = 3
    ENS_ID = 4
    GENE_NAME = 5
    GENE_TYPE = 6

class GenCode:
    """-gencode- information. Yes, without 'e' in 'gene'"""
    def __init__(self):
        self.ensemble_version = None
        self.ensemble = None
        self.version = None
        self.name = None
        self.start_location = None
        self.end_location = None

    @classmethod
    def loadFromShortRow(cls, row):
        gencode = GenCode()
        gencode.ensemble_version = row[GFTFS.ENS_ID]
        gencode.ensemble = gencode.ensemble_version.split('.')[0]
        gencode.version = gencode.ensemble_version.split('.')[1]
        gencode.name = row[GFTFS.GENE_NAME]

    @classmethod
    def loadFromGTFRow(cls, comps):
        gencode = GenCode()

        key = None
        value = None
        key_value_pairs = [x.translate(None,';') for x in comps[GFTF.KEY_VALUE_PAIRS:]]

        for i,string in enumerate(key_value_pairs):
            if key is None:
                key = string
            elif value is None:
                value = string.translate(None,'"\n')
                if key == GFTF.GENE_ID:
                    gencode.ensemble_version = value.translate(None,'"')
                    gencode.ensemble = gencode.ensemble_version.split('.')[0]
                    gencode.version = gencode.ensemble_version.split('.')[1]
                elif key == GFTF.GENE_NAME:
                    gencode.name = value
                key = None
                value = None

        return gencode

def parse_gencode_file(path, callback):
    class Wrapper(object):
        def __init__(self, c):
            self.callback = c

        def __call__(self, i, comps):
            if "##" in comps[0]:
                return

            feature = comps[GFTF.FEATURE_TYPE]
            if feature in K_NOT_GENES:
                return

            gencode = GenCode.loadFromGTFRow(comps)
            self.callback(gencode)

    wrapper = Wrapper(callback)
    Utilities.parse_file(path, wrapper)

class GenCodeSet:
    """Collection of gencode data sets"""
    def __init__(self):
        self.gencodes = []
        self.gencodes_by_ensemble_id = {}
        self.gencodes_by_ensemble_version = {}

    @classmethod
    def LoadGeneCodeInput(cls, file_name):
        gencodes = GenCodeSet()
        with open(file_name, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
            for row in reader:
                gencode = GenCode.loadFromShortRow(row)
                gencodes.gencodes.append(gencode)

                if not gencode.ensemble in gencodes.gencodes_by_ensemble_id:
                    gencodes.gencodes_by_ensemble_id[gencode.ensemble] = gencode
                else:
                    raise Exception('Duplicate ensemble id, check'+gencode.ensemble)

                if not gencode.ensemble_version in gencodes.gencodes_by_ensemble_version:
                    gencodes.gencodes_by_ensemble_version[gencode.ensemble_version] = gencode
                else:
                    raise Exception('Duplicate ensemble version, check'+gencode.ensemble_version)
        return gencodes

    @classmethod
    def LoadGTF(cls, file_name):
        gencodes = GenCodeSet()
        with open(file_name, 'rb') as tabfile:
            for line in tabfile:
                if "##" in line:
                    continue

                row = line.split("\t")
                #print row
                feature = row[GFTF.FEATURE_TYPE]
                if feature in K_NOT_GENES:
                    continue

                gencode = GenCode.loadFromGTFRow(row)

                gencodes.gencodes.append(gencode)

                if gencode.ensemble in gencodes.gencodes_by_ensemble_id:
                    raise Exception('Duplicate ensemble id, check'+gencode.ensemble)
                if gencode.ensemble_version in gencodes.gencodes_by_ensemble_version:
                    raise Exception('Duplicate ensemble version, check'+gencode.ensemble_version)

                gencodes.gencodes_by_ensemble_id[gencode.ensemble] = gencode
                gencodes.gencodes_by_ensemble_version[gencode.ensemble_version] = gencode

        return gencodes