import gzip
import sqlite3
import logging
import os
import math

def parse_file(path, callback, expect_throw=False):
    def parse(path, callback):
        if ".gz" in path:
            with gzip.open(path) as file:
                for i, line in enumerate(file):
                    callback(i, line.strip().split())
        else:
            with open(path) as file:
                for i, line in enumerate(file):
                    callback(i, line.strip().split())


    if expect_throw:
        try:
            parse(path, callback)
        except IOError:
            logging.info("Expected bad file, and that's what we got: %s" % (path, ))
        except Exception as e:
            raise e
    else:
        parse(path, callback)


def connect(db_path):
    if os.path.exists(db_path):
        raise RuntimeError("Existing db")
    return sqlite3.connect(db_path)

def setup_db(connection):
    cursor = connection.cursor()
    cursor.execute("CREATE TABLE extra(gene TEXT, genename TEXT, `n.snps.in.model` INTEGER, `pred.perf.R2` DOUBLE, `pred.perf.pval` DOUBLE, `pred.perf.qval`DOUBLE);")
    cursor.execute("CREATE TABLE weights(rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER);")
    cursor.execute("CREATE INDEX extra_gene ON extra(gene);")
    cursor.execute("CREATE INDEX weights_gene ON weights(gene);")
    cursor.execute("CREATE INDEX weights_rsid ON weights(rsid);")
    cursor.execute("CREATE INDEX weights_rsid_gene ON weights(rsid, gene);")
    connection.commit()

# Weight DB input row format. So that we know how to read in data.
class WDBIF(object):
    SNP = 0
    GENE = 1
    GENE_NAME = 2
    REFERENCE_ALLELE = 3
    EFFECT_ALLELE = 4
    WEIGHT = 5
    N_SNP = 6
    GENE_R2 = 7
    GENE_PVALUE = 8
    GENE_QVALUE = 9

# Send tuples because they use less memory. Fullfledged objects might be too much.
def insert_entries(db, gene_entries):
    cursor = db.cursor()
    keys = sorted(gene_entries.keys())
    i = 0
    for key in keys:
        rows = gene_entries[key]
        if len(rows) == 0:
            continue
        i += len(rows)
        insert = [(e[WDBIF.SNP], e[WDBIF.GENE], e[WDBIF.WEIGHT], e[WDBIF.REFERENCE_ALLELE], e[WDBIF.EFFECT_ALLELE],) for e in rows]
        cursor.executemany("INSERT INTO weights(rsid, gene, weight, ref_allele, eff_allele) VALUES(?, ?, ?, ?, ?)", insert)
    logging.info("Inserted %d snp entries", i)

    logging.info("Inserting gene entries")
    i = 0
    for gene, rows in gene_entries.iteritems():
        if len(rows) == 0:
            continue
        r = rows[0]
        i += 1
        insert = (r[WDBIF.GENE], r[WDBIF.GENE_NAME], len(rows), r[WDBIF.GENE_R2], r[WDBIF.GENE_PVALUE], r[WDBIF.GENE_QVALUE])
        cursor.execute("INSERT INTO extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) VALUES(?, ?, ?, ?, ?, ?)", insert)
    logging.info("Inserted %d gene entries", i)
    db.commit()

import gzip
class FileIterator(object):
    def __init__(self, path, header=None, compressed = False, ignore_until_header = False):
        self.path = path
        self.compressed = compressed
        self.header = header
        self.ignore_until_header = ignore_until_header
        if ignore_until_header and not header:
            raise RuntimeError("File iterator received conflicting header information")

    def iterate(self,callback=None):
        if self.compressed:
            with gzip.open(self.path, 'rb') as file_object:
                self._iterateOverFile(file_object, callback)
        else:
            with open(self.path, 'rb') as file_object:
                self._iterateOverFile(file_object, callback)

    def _iterateOverFile(self, file_object, callback):
        if self.ignore_until_header:
            self._ignore_until_header(file_object)
        else:
            if self.header is not None:
                line = file_object.readline().strip("\n")
                if len(self.header) and line != self.header:
                    raise RuntimeError(self.path, "Unexpected header: got:\n%s\n,expected:%s\n" %(line, self.header))

        self._processFile(file_object, callback)

    def _ignore_until_header(self, file_object):
        if self.ignore_until_header and self.header:
            skip = True
            while skip:
                l = file_object.readline()
                if not l:
                    raise RuntimeError("Wrong header lookup in %s" % (self.path,))
                l = l.strip()
                if self.header in l:
                    skip = False


    def _processFile(self, file_object, callback):
        if callback is not None:
            for i,line in enumerate(file_object):
                callback(i, line)

import csv
class CSVFileIterator(FileIterator):
    def __init__(self, path, header=None, compressed = False, ignore_until_header = False, delimiter = " ", quotechar='"'):
        super(CSVFileIterator, self).__init__(path, header, compressed, ignore_until_header)
        self.delimiter = delimiter
        self.quotechar = quotechar

    def _processFile(self, file_object, callback):
        if callback is not None:
            reader = csv.reader(file_object, delimiter=self.delimiter, quotechar=self.quotechar)
            for i,row in enumerate(reader):
                callback(i, row)

#keep all snps
def process_all_rows(row, genes, gene_key):
    gene_key = row[gene_key]
    if not gene_key in genes:
        genes[gene_key] = {}
    rows = genes[gene_key]
    snp = row[WDBIF.SNP]
    if snp in rows:
        #overwrite existing only if better
        existing = rows[snp]
        if math.fabs(float(existing[WDBIF.WEIGHT])) < math.fabs(float(row[WDBIF.WEIGHT])):
            rows[snp] = row
    else:
        rows[snp] = row

#keep only the best snp weight
def keep_best_row(row, genes, gene_key):
    gene_key = row[gene_key]
    if not gene_key in genes:
        genes[gene_key] = {}

    rows = genes[gene_key]
    snp = row[WDBIF.SNP]
    if len(rows) == 0:
        rows[snp] = row
    else:
        existing_key = rows.keys()[0]
        r = rows[existing_key]
        if math.fabs(float(r[WDBIF.WEIGHT])) < math.fabs(float(row[WDBIF.WEIGHT])):
            del rows[existing_key]
            rows[snp] = row
