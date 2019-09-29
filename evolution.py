import argparse
import csv

import numpy as np

parser = argparse.ArgumentParser(description="Reload a results file and track the minimum.")

parser.add_argument('source',   help='Specify one or more input file names.',     type=str, nargs='+')
parser.add_argument('--out',    help='Specify output file name [evolution.csv].', type=str, default='evolution.csv')
parser.add_argument('--column', help='Specify column of interest [4 (cost)].',    type=int, default=4)
parser.add_argument('--rows',   help='Number of data points to read [1000].',     type=int, default=1000)

args = parser.parse_args()

column = args.column - 1
Nfiles = len(args.source)
result = np.zeros((args.rows, 4+Nfiles))

cid = 4
for file_name in args.source:
    print('Source file: {s}'.format(s=file_name))
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rid = 0
        count = 0
        for row in reader:
            rid = rid + 1
            if rid > 2:
                values = np.asarray(list(map(float, row)))

                if count == 0:
                    result[count,cid] = values[column]
                elif values[column] < result[count-1,cid]:
                    result[count,cid] = values[column]
                else:
                    result[count,cid] = result[count-1,cid]

                count += 1
                if count >= args.rows:
                    break
    cid += 1

result[:,0] = np.average(result[:,4:], axis=1)
result[:,1] = np.amin(result[:,4:], axis=1)
result[:,2] = np.amax(result[:,4:], axis=1)
result[:,3] = np.std(result[:,4:], axis=1)

with open(args.out, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['Mean','Minimum','Maximum','Std. Dev.'])
    for i in range(0, args.rows):
        writer.writerow(result[i,0:4])
