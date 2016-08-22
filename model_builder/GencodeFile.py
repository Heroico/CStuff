__author__ = 'heroico'
#trimmed from PredictDBAnalysis/gencode_input

import csv
import gzip
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
    GENOMIC_STRAND = 6
    GENOMIC_PHASE = 7
    KEY_VALUE_PAIRS = 8

    #there are several other key-value pairs but we are concerned with these
    GENE_ID = "gene_id"
    TRANSCRIPT_ID = "transcript_id"
    GENE_TYPE = "gene_type"
    GENE_STATUS = "gene_status"
    GENE_NAME = "gene_status"
    TRANSCRIPT_TYPE = "transcript_type"
    TRANSCRIPT_STATUS = "transcript_status"
    TRANSCRIPT_NAME = "transcript_name"
    EXON_NUMBER = "exon_number"
    EXON_ID = "exon_id"
    LEVEL = "level"

    #some are missing
    TAG = "tag"


class GFTFS:
    """gencode short file table format"""
    CHROMOSOME = 0
    FEATURE_TYPE = 1
    N_START_LOCATION = 2
    N_END_LOCATION = 3
    ENS_ID = 4
    GENE_NAME = 5
    GENE_TYPE = 6

class GenCode(object):
    """-gencode- information. Yes, without 'e' in 'gene'"""
    def __init__(self):
        #fixed fields
        self.chromosome_name = None
        self.annotation_source = None
        self.feature_type = None
        self.start_location = None
        self.end_location = None
        self.score = None
        self.genomic_strand = None
        self.genomic_phase = None

        #key-value pairs, mandatory
        self.gene_id = None
        self.transcript_id = None
        self.gene_type = None
        self.gene_status = None
        self.gene_name = None
        self.transcript_type = None
        self.transcript_status = None
        self.transcript_name = None
        self.exon_number = None
        self.exon_id = None
        self.level = None

        #key-value pairs, mandatory
        self.tag = None
        self.ccdsid = None
        self.havana_gene = None
        self.havana_transcript = None
        self.protein_id = None
        self.ont = None
        self.transcription_support_level = None

        self.ensemble_version = None
        self.ensemble = None
        self.version = None
        self.name = None

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

        gencode.chromosome_name = comps[GFTF.CHROMOSOME]
        gencode.annotation_source = comps[GFTF.ANNOTATION_SOURCE]
        gencode.feature_type = comps[GFTF.FEATURE_TYPE]
        gencode.start_location = comps[GFTF.N_START_LOCATION]
        gencode.end_location = comps[GFTF.N_END_LOCATION]
        gencode.score = comps[GFTF.SCORE]
        gencode.genomic_strand = comps[GFTF.GENOMIC_STRAND]
        gencode.genomic_phase = comps[GFTF.GENOMIC_PHASE]

        key = None
        value = None
        key_value_pairs = [x.translate(None,';') for x in comps[GFTF.KEY_VALUE_PAIRS:]]

        for i,string in enumerate(key_value_pairs):
            if key is None:
                key = string
            elif value is None:
                value = string.translate(None,'"\n')
                if key == GFTF.GENE_ID:
                    gencode.gene_id = value
                    gencode.ensemble_version = value.translate(None,'"')
                    gencode.ensemble = gencode.ensemble_version.split('.')[0]
                    gencode.version = gencode.ensemble_version.split('.')[1]
                elif key == GFTF.TRANSCRIPT_ID:
                    gencode.transcript_id = value
                elif key == GFTF.GENE_TYPE:
                    gencode.gene_type = value
                elif key == GFTF.GENE_STATUS:
                    gencode.gene_status = value
                elif key == GFTF.GENE_NAME:
                    gencode.gene_name = value
                    gencode.name = value
                elif key == GFTF.TRANSCRIPT_TYPE:
                    gencode.transcript_type = value
                elif key == GFTF.TRANSCRIPT_STATUS:
                    gencode.transcript_status = value
                elif key == GFTF.TRANSCRIPT_NAME:
                    gencode.transcript_name = value
                elif key == GFTF.EXON_NUMBER:
                    gencode.exon_number = value
                elif key == GFTF.EXON_ID:
                    gencode.exon_id = value
                elif key == GFTF.LEVEL:
                    gencode.exon_id = value
                elif key == GFTF.TAG:
                    gencode.tag = value
                key = None
                value = None

        return gencode


class GathererCallback(object):
    def __init__(self):
        self.results = []
    def __call__(self, gencode):
        self.results.append(gencode)

def parse_gencode_file(path, callback, only_genes=True):
    class Wrapper(object):
        def __init__(self, c):
            self.callback = c

        def __call__(self, i, comps):
            if "##" in comps[0]:
                return

            feature = comps[GFTF.FEATURE_TYPE]
            if only_genes and feature in K_NOT_GENES:
                return

            gencode = GenCode.loadFromGTFRow(comps)
            self.callback(gencode)

    wrapper = Wrapper(callback)
    Utilities.parse_file(path, wrapper)

class BuildGeneNameRelationShipCallback(object):
    def __init__(self):
        self.ensemble_to_name = {}
        self.name_to_ensemble = {}

    def __call__(self, i, comps):
        if "##" in comps[0]:
            return

        feature = comps[GFTF.FEATURE_TYPE]
        if feature != "gene":
            return

        gene_id = None
        gene_name =None
        key = None
        value = None
        key_value_pairs = [x.translate(None,';') for x in comps[GFTF.KEY_VALUE_PAIRS:]]

        for i,string in enumerate(key_value_pairs):
            if key is None:
                key = string
            elif value is None:
                value = string.translate(None,'"\n')
                if key == GFTF.GENE_ID:
                    gene_id = value.translate(None,'"')
                elif key == GFTF.GENE_NAME:
                    gene_name = value.translate(None,'"')
                key = None
                value = None

        self.ensemble_to_name[gene_id] = gene_name
        self.name_to_ensemble[gene_name] = gene_id

def ensemble_to_name_relationships(path):
    callback = BuildGeneNameRelationShipCallback()
    Utilities.parse_file(path, callback)
    return callback.ensemble_to_name, callback.name_to_ensemble


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