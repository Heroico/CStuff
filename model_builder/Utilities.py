import gzip
import logging

def parse_file(path, callback):
    with gzip.open(path) as file:
        try:
            for i, line in enumerate(file):
                callback(i, line.strip().split())
        except IOError:
            logging.info("Expected bad file, and that's what we got: %s" % (path, ))
        except Exception as e:
            raise e