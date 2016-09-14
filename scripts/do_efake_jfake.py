#! /usr/bin/env python

import os
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description='do_efake_jfake')

parser.add_argument('-v', dest='version', required=True, help='Mini version')
parser.add_argument('-d', dest='dids', help='Specific run (comma separated)')

args = parser.parse_args()

mini_dir  = '/raid/falonso/mini2/'

data15 = glob.glob(mini_dir+'v%s/data15_13TeV.*' % args.version) 
data16 = glob.glob(mini_dir+'v%s/data16_13TeV.*' % args.version) 

all_data = data15 + data16

if args.dids is not None:

    dids = args.dids.split(',')
    
    only_data = []
    for s in all_data:
        for did in dids:
            if did in s:
                only_data.append(s)

    all_data = only_data


for s in sorted(all_data):

    print 'Processing', s

    output_path_efake = s.replace('data15_13TeV', 'efake15').replace('data16_13TeV', 'efake16')
    output_path_jfake = s.replace('data15_13TeV', 'jfake15').replace('data16_13TeV', 'jfake16')

    # efake
    cmd = 'create_efake_mini %s %s' % (s, output_path_efake)
    os.system(cmd)

    # jfake
    cmd = 'create_jfake_mini %s %s' % (s, output_path_jfake)
    os.system(cmd)


