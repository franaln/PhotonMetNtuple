#! /usr/bin/env python

import os
import sys
from samples import all_samples

#xs_file = os.environ['ROOTCOREDIR'] + '/data/SUSYTools/susy_crosssections_13TeV.txt'
xs_file = sys.argv[1]

lines = open(xs_file).read().split('\n')

def get_xs(did):
    
    
    for line in lines:
        if not line or line.startswith('#'):
            continue
    
        if line.startswith(did):
    
            did, name, xs, kfact, eff, relunc = line.split()
    
            return float(xs)*float(kfact)*float(eff)
    
    print 'xs not found for %s' % did



with open('xs_dict.py', 'w+') as of:

    of.write('xs_dict = dict()\n\n')

    for sample in all_samples:

        if 'data' in sample:
            continue

        did = sample.split('.')[1]

        xs = get_xs(did)

        print did, xs
        of.write('xs_dict[%s] = %s\n' % (did, xs))

