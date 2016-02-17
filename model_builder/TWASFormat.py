import numpy
import os

def build_subpaths(folder, subfolder):
    sub_path = os.path.join(folder, subfolder)
    sub_path = os.path.join(sub_path, subfolder)
    return sub_path

class MTF(object):
    snp=0
    pos=1
    a1=2
    a2=3

def load_map(path):
    snps = []
    with open(path) as file:
        for i, line in enumerate(file):
            if i==0:
                continue
            comps = line.strip().split()
            row = (comps[0], comps[1], comps[2], comps[3], )
            snps.append(row)
    return snps

def load_cor(path):
    cors = []
    with open(path) as file:
        for line in file:
            c = line.strip()
            cors.append(float(c))
    return cors

def load_ld(path):
    rows = []
    with open(path) as file:
        for line in file:
            row = [float(x) for x in line.strip().split()]
            rows.append(row)
    array = numpy.array(rows)
    i, j = numpy.indices(array.shape)
    array[i == j] += 0.01
    return array

def build_weights(sub_path):
    cor_path = sub_path + ".wgt.cor"
    cors = load_cor(cor_path)

    ld_path = sub_path + ".wgt.ld"
    ld = load_ld(ld_path)

    weights = calculate_weights(cors, ld)
    return weights

def calculate_weights(cors, ld):
    inv = numpy.linalg.inv(ld)
    dot = numpy.dot(cors, inv)
    return dot