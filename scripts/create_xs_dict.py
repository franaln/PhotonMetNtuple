#! /usr/bin/env python

import os
import sys

#xs_file = os.environ['ROOTCOREDIR'] + '/data/SUSYTools/susy_crosssections_13TeV.txt'
try:
    xs_file = sys.argv[1]
    sample_files = [ i for i in sys.argv[2:] ]
except:
    print 'usage: create_xs_dict.py [xs-file] [sample-file] [sample-file2] ...'
    sys.exit(1)


xs_lines = open(xs_file).read().split('\n')

slines = []
for sfile in sample_files:
    slines.extend(open(sfile).read().split('\n'))


def get_xs(did):
    
    for line in xs_lines:
        if not line or line.startswith('#'):
            continue
    
        if line.startswith(did):
    
            did, name, xs, kfact, eff, relunc = line.split()
    
            return float(xs)*float(kfact)*float(eff)
    
    print 'xs not found for %s' % did



xs_dict = dict()

for line in slines:

    if not line or line.startswith('#'):
        continue

    did = line.split('.')[1]

    xs_dict[did] = get_xs(did)



for did, xs in sorted(xs_dict.iteritems()):

    print 'xs_dict[%s] = %s' % (did, xs)

