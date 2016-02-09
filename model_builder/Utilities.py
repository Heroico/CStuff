import gzip
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
