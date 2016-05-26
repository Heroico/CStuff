import gzip
import sqlite3
import logging
import os

def parse_file(path, callback, expect_throw=False):
    def parse(path, callback):
        with gzip.open(path) as file:
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
    cursor.execute("CREATE TABLE extra (gene TEXT, genename TEXT, R2 DOUBLE,  `n.snps` INTEGER, pval DOUBLE)")
    cursor.execute("CREATE INDEX extra_gene ON extra (gene)")
    cursor.execute("CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
    cursor.execute("CREATE INDEX weights_rsid ON weights (rsid)")
    cursor.execute("CREATE INDEX weights_gene ON weights (gene)")
    cursor.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
    connection.commit()

class WDBIF(object):
    SNP = 0
    GENE = 1
    GENE_NAME = 2
    WEIGHT = 3
    REFERENCE_ALLELE = 4
    EFFECT_ALLELE = 5

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
        cursor.executemany("INSERT INTO weights VALUES(?, ?, ?, ?, ?, NULL, NULL, NULL)", insert)
    logging.info("Inserted %d snp entries", i)

    logging.info("Inserting gene entries")
    i = 0
    for gene, rows in gene_entries.iteritems():
        if len(rows) == 0:
            continue
        r = rows[0]
        i += 1
        insert = (e[WDBIF.GENE], e[WDBIF.GENE_NAME], "NA", len(rows), "NA")
        cursor.execute("INSERT INTO extra VALUES(?, ?, ?, ?, ?)", insert)
    logging.info("Inserted %d gene entries", i)
    db.commit()
