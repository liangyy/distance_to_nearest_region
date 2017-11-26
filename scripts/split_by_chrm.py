import argparse
parser = argparse.ArgumentParser(prog='split_by_chrm.py')
parser.add_argument('--chrm')
parser.add_argument('--input')
parser.add_argument('--output')


args = parser.parse_args()

import os
chrm = args.chrm
cmd = 'zcat {inp} | grep {chrm} | gzip > {output}'.format(inp = args.input,
    chrm = chrm,
    output = args.output)
os.system(cmd)
