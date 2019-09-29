import argparse
import csv

import numpy as np

from BeesEtAl.Base_Sorter import Base_Sorter

parser = argparse.ArgumentParser(description="Reload a results file and determine the Pareto-optimal set.")

parser.add_argument('source',          help='Specify one or more input file names.',                    type=str, nargs='+')
parser.add_argument('--out',           help='Specify output file name [pareto.csv].',                   type=str, default='pareto.csv')
parser.add_argument('--pareto-in',     help='The input files are Pareto files.',                                                    action='store_true')
parser.add_argument('--without-Nc',    help='The input files are old-style results files without number of circumferential rings.', action='store_true')
parser.add_argument('--filter-hybrid', help='Filter input results for hybrid designs (for new-style results files only).',          action='store_true')
parser.add_argument('--filter-rings',  help='Filter input results for specified number of node-rings.', type=int, default=0, dest='rings')
parser.add_argument('--row-limit',     help='Only read the specified number of rows.',                  type=int, default=0)

args = parser.parse_args()

S = None

for file_name in args.source:
    print('Source file: {s}'.format(s=file_name))
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rid = 0
        count = 0
        for row in reader:
            rid = rid + 1
            if rid > 2 or (rid > 1 and args.pareto_in):
                values = np.asarray(list(map(float, row)))
                if args.pareto_in:
                    if len(values) == 30: # oops - missing the Nr column
                        X    = np.copy(values[1:])
                        X[0] = args.rings
                    else:
                        X    = values[2:]
                    cost   = values[0:2]
                elif args.without_Nc:
                    X      = values[4:]
                    cost   = np.asarray([values[3], values[2] * 1E3])
                else:
                    if args.filter_hybrid and values[5] - values[4] < 1.5:
                        # For purely lamella designs, Nc = Nr - 1
                        continue
                    X      = values[5:]
                    cost   = np.asarray([values[3], values[2] * 1E3])

                if args.rings > 0:
                    if int(X[0]) != args.rings:
                        continue

                if S is None:
                    S = Base_Sorter(len(X))
                S.push(cost, X)

                count = count + 1
                if args.row_limit > 0 and count >= args.row_limit:
                    break

S.pareto(args.out)
