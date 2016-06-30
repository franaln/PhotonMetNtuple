#! /usr/bin/env python

import os
import sys

#xs_file = os.environ['ROOTCOREDIR'] + '/data/SUSYTools/susy_crosssections_13TeV.txt'
try:
    xs_file = sys.argv[1]
    sample_file = sys.argv[2]
except:
    print 'usage: create_xs_dict.py [xs_file] [sample_file]'
    sys.exit(1)


xs_lines = open(xs_file).read().split('\n')
sample_lines = open(sample_file).read().split('\n')

def get_xs(did):
    
    for line in xs_lines:
        if not line or line.startswith('#'):
            continue
    
        if line.startswith(did):
    
            did, name, xs, kfact, eff, relunc = line.split()
    
            return float(xs)*float(kfact)*float(eff)
    
    print 'xs not found for %s' % did



with open('xs_dict.py', 'w+') as of:

    of.write('xs_dict = dict()\n\n')

    for line in sample_lines:

        if not line or line.startswith('#'):
            continue

        did = line.split('.')[1]

        xs = get_xs(did)

        print 'xs for %s: %s' % (did, xs)
        of.write('xs_dict[%s] = %s\n' % (did, xs))

