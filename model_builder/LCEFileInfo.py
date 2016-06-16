import gzip
import logging

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
