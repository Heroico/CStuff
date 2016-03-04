import gzip

def load_variance(variance_path):
    vars = {}
    with gzip.open(variance_path, "rb") as var_file:
        for line in var_file:
            comps = line.strip().split(",")
            vars[comps[0]] = float(comps[1])
    return vars