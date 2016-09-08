#! /usr/bin/env python
import logging
import sys
import os
import csv
import re
regexp = re.compile('[_]')

def fix_row(row, last_row,  args):
    def f(x):
        x = regexp.sub(" ", x)

        if args.decimals and "." in x:
            try:
                number = float(x)
                x = "%.2E" % number
                raise RuntimeError("success")
            except:
                pass

        return x
    row = [f(x) for x in row]

    up_to = args.clear_repeated_column_up_to
    if up_to and last_row:
        if row[0: up_to] == last_row[0:up_to]:
            row = ["" if i<up_to else x for i,x in enumerate(row)]

    return row

def run(args):
    if os.path.exists(args.output_path):
        logging.info("Output already there, delete it if you want it done again")
        return

    folder = os.path.split(args.output_path)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open(args.input_path) as in_file:
        with open(args.output_path, "w") as out_file:
            reader = csv.reader(in_file, delimiter = args.delimiter)
            last_row = None
            for i, row in enumerate(reader):
                fixed_row = fix_row(row, last_row, args)
                last_row = row
                if i == 0:
                    header = "\\begin{{tabular}}{{ {} }}\n".format( " | ".join(["c" for  x in fixed_row]))
                    out_file.write(header)
                line = "{} \\\\\n".format(" & ".join(fixed_row))
                out_file.write(line)
            out_file.write("\\end{tabular}")


def configureLogging(level=5, target=sys.stdout):
    logger = logging.getLogger()
    logger.setLevel(level)

    # create console handler and set level to info
    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert CSV ot TSV to table format")

    parser.add_argument("--input_path",
                        help="path to delimited file",
                        default=None)

    parser.add_argument("--output_path",
                        help="path to output",
                        default=None)

    parser.add_argument("--delimiter",
                        help= "snp metadata file",
                        default=None)

    parser.add_argument("--decimals",
                        help= "round any numbers to this many decimals",
                        type=int,
                        default=None)

    parser.add_argument("--clear_repeated_column_up_to",
                        help= "Empty components of a row if they match the last one's.",
                        type=int,
                        default=None)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    configureLogging(int(args.verbosity))

    run(args)