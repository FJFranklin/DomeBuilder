import argparse
import csv

import numpy as np

from BeesEtAl.Base_Sorter import Base_Sorter

parser = argparse.ArgumentParser(description="Reload a results file and determine the Pareto-optimal set.")

parser.add_argument('source',         help='Specify one or more input file names.', type=str, nargs='+')
parser.add_argument('--out',          help='Specify output file name [pareto.csv].', type=str, default='pareto.csv')
parser.add_argument('--pareto-in',    help='The input files are Pareto files.', action='store_true')
parser.add_argument('--without-Nc',   help='The input files are old-style results files without number of circumferential rings.', action='store_true')
parser.add_argument('--filter-rings', help='Filter input results for specified number of node-rings.', dest='rings', default=0, type=int)

args = parser.parse_args()

S = None

for file_name in args.source:
    print('Source file: {s}'.format(s=file_name))
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rid = 0
        for row in reader:
            rid = rid + 1
            if rid > 2 or (rid > 1 and args.pareto_in):
                values = np.asarray(list(map(float, row)))
                if args.pareto_in:
                    X      = values[2:]
                    cost   = values[0:2]
                elif args.without_Nc:
                    X      = values[4:]
                    cost   = np.asarray([values[3], values[2] * 1E3])
                else:
                    X      = values[5:]
                    cost   = np.asarray([values[3], values[2] * 1E3])

                if args.rings > 0:
                    if int(X[0]) != args.rings:
                        continue

                if S is None:
                    S = Base_Sorter(len(X))
                S.push(cost, X)

S.pareto(args.out)
