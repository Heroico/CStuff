import gzip
import sqlite3
import logging

def parse_file(path, callback, trow=False):
    def parse(path, callback):
        with gzip.open(path) as file:
            for i, line in enumerate(file):
                callback(i, line.strip().split())

    if trow:
        try:
            parse(path, callback)
        except IOError:
            logging.info("Expected bad file, and that's what we got: %s" % (path, ))
        except Exception as e:
            raise e
    else:
        parse(path, callback)

def connect(db_path):
    return sqlite3.connect(db_path)

def setup_db(connection):
    cursor = connection.cursor()
    cursor.execute("CREATE TABLE extra (gene TEXT, genename TEXT, R2 DOUBLE,  `n.snps` INTEGER)")
    cursor.execute("CREATE INDEX extra_gene ON extra (gene)")
    cursor.execute("CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
    cursor.execute("CREATE INDEX weights_rsid ON weights (rsid)")
    cursor.execute("CREATE INDEX weights_gene ON weights (gene)")
    cursor.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
    connection.commit()

def insert_entries(db, genes):
    cursor = db.cursor()
    keys = sorted(genes.keys())
    i = 0
    for key in keys:
        rows = genes[key]
        if len(rows) == 0:
            continue
        i += len(rows)
        cursor.executemany("INSERT INTO weights VALUES(?, ?, ?, ?, ?, NULL, NULL, NULL)", [(r[0], r[1], r[3], r[4], r[5]) for r in rows])
    logging.info("Inserted %d snp entries", i)

    logging.info("Inserting gene entries")
    i = 0
    for gene, rows in genes.iteritems():
        if len(rows) == 0:
            continue
        r = rows[0]
        i += 1
        cursor.execute("INSERT INTO extra VALUES(?, ?, ?, ?)", (r[1], r[2], "NA", len(rows)))
    logging.info("Inserted %d gene entries", i)
    db.commit()
