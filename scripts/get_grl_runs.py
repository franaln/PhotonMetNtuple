#! /usr/bin/env python

import sys
from xml.etree.ElementTree import ElementTree

if len(sys.argv) < 2:
    print 'usage: get_grl_runs.pt [grl-file]'
    sys.exit()

tree = ElementTree()
tree.parse(sys.argv[1]) 

links = list(tree.iter("Metadata")) 

grl_runs = []
for l in links:
    if l.get('Name') == 'RunList':
        grl_runs = l.text.split(',')


for run in grl_runs:
    print run
