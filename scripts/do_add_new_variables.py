#! /usr/bin/env python

import os
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description='do_efake_jfake')

parser.add_argument('--vin',  required=True, help='Input version')
parser.add_argument('--vout', required=True, help='Output version')
parser.add_argument('--filter', help='Filter samples')

args = parser.parse_args()

mini_dir  = '/raid/falonso/mini2/'

samples = [ p for p in glob.glob(mini_dir+'v%s/*.root' % args.vin) if os.path.isfile(p) ]

if args.filter is not None:
    samples = [ s for s in samples if args.filter in s ]

for s in sorted(samples):

    print 'Processing', s

    output_path = s.replace('v%s' % args.vin, 'v%s' % args.vout)

    cmd = 'add_new_variables %s %s' % (s, output_path)
    os.system(cmd)

