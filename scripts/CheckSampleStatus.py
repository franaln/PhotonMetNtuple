#!/usr/bin/env python

##############################################################
#                                                            #
# CheckSampleStatus.py                                       #
# C. Ohm, M. Tripiana, 2016                                  #
#                                                            #
# Script for finding out the status of a list of samples     #
# defined by their DSIDs and a list of AMI tags              #
#                                                            #
##############################################################

# Set the default tags you want to use for your sample status check
# For MC15c, see https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AtlasProductionGroupMC15c
defaultRecoTagFilteringPattern = "r7" # this is for the rucio query, which is filtered by the below
defaultSimTags = ["s2608", "s2726"] # MC15 standard tags
defaultRecoTags = ["f*", "r*", "r7725", "a818", "a821"] # for AF-II, enable this line
defaultDerivationTags = ["p2769", "p2709", "p2689","p2666", "p2667", "p2645", "p2622", "p2623", "p2613", "p2614", ] # in order of priority (in case several are available)

import sys, os, argparse, subprocess, re, pickle
from collections import OrderedDict

def runCommand(cmd, verbose = False):
    if verbose:
        print "  Will run the following command: %s" % cmd
    cmdResult = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return cmdResult.stdout

def getSamplesFromPattern(dsid,
                         projectTag = "mc15_13TeV",
                         simTags = defaultSimTags,
                         recoTags = defaultRecoTags,
                         verbose = False,
                         format = "DAOD_SUSY1",
                         pTags = [],
                         derivation = ""):
    
    pattern = projectTag + '.' + dsid +'.physics_Main*AOD*'

    # do the rucio query
    queryPattern = pattern
    cmd = "rucio ls --short --filter type=CONTAINER %s | sort -r " % (queryPattern) # sort -r prioritizes higher e-tags
    queryResult = runCommand(cmd, verbose)

    # preprocess the output
    lines = [line.rstrip().replace(projectTag+":", "") for line in queryResult.readlines()]

    dsPattern = pattern.replace(defaultRecoTagFilteringPattern+"*", "").replace("*", ".*")
    dsRE = re.compile(dsPattern)
    datasets = filter(dsRE.match, lines)
    if verbose:
        print "  Dataset regex: %s" % dsPattern
        print "  All datasets: "
        for ds in datasets:
            print "    %s" % ds

    aodPattern = dsPattern.replace("AOD.*", "\.merge\.AOD\..*")
    aodPattern += "("+'|'.join(simTags)+")?.*("+'|'.join(recoTags)+")?.*"
    aodRE = re.compile(aodPattern)
    aodDatasets = filter(aodRE.match, datasets)
    if verbose:
        print "  AOD regex: %s" % aodPattern
        print "  AOD datasets: "
        for ds in aodDatasets:
            print "    %s" % ds

    if not aodDatasets:
        aodDatasets.append(" N/A ")

    daodPattern = aodPattern.replace("AOD", format)
    daodPattern += "("+'|'.join(pTags)+")"
    daodRE = re.compile(daodPattern)
    daodDatasets = filter(daodRE.match, datasets)
    if verbose:
        print "  DAOD regex: %s" % daodPattern
        print "  DAOD datasets:"
        for ds in daodDatasets:
            print "    %s" % ds
    # if there are several matching the requested tags, pick the one with the tag mentioned first in the list
    if len(daodDatasets) > 1:
        print "Found more than one:"
        for ds in daodDatasets:
            print "  %s" % ds
        for pTag in pTags:
            if any(ds for ds in daodDatasets if pTag in ds):
                daodDatasets = [next(ds for ds in daodDatasets if pTag in ds)]
                break # we're done after we've found the first
        print "Was the right one selected?"
        print "  %s" % daodDatasets[0]
    
    # if there are no matching DAOD datasets, fill a dummy string
    if len(daodDatasets) < 1:
        daodDatasets.append(" N/A ")

    # check for non-dataset line in the output - they could be e.g. error messages
    otherLines = [x for x in lines if x not in datasets] 
    if len(otherLines) > 0:
        print "Other lines were found in the query output - there may have been a problem:"
        for line in otherLines:
            print line

    returnDict = {"AOD": aodDatasets[0], "DAOD": daodDatasets[0]}
    return returnDict

def getAmiStatus(ds, verbose):
    cmd = "ami show dataset info %s | grep prodsysStatus" % ds
    amiResult = runCommand(cmd, verbose)


    # try:
    #     cmd = 'ami show dataset info %s | grep totalEvents' % ds
    #     totalEvents = runCommand(cmd, verbose).readlines()[0].split(':')[1].strip()
    # except:
    #     totalEvents = '0'

    lines = [line.rstrip() for line in amiResult.readlines()]
    if verbose:
        for line in lines:
            print line
    if len(lines) == 1:
        amiStatus = lines[0].split(':')[1].strip()
    elif len(lines) > 1:
        print "Weird AMI status for %s (saving N/A):" % ds
        print lines
        amiStatus = "N/A"
    else:
        print "No AMI status for %s" % ds
        amiStatus = "N/A"

    #amiStatus = '%s (%s)' % (amiStatus, totalEvents)

    return amiStatus

def printSamplesDict(s):
    row_format ="{0:<15}{1:^30}{2:^30}{3:<120}"
    print row_format.format("=== DSID ===", "=== AOD ===", "=== DAOD (AOD) ===", "=== DAOD dataset name ===")
    for ds in s:
        dsName = s[ds]["DAOD"]
        if dsName == " N/A ":
            dsName = "(AOD: " + s[ds]["AOD"] + ")"
        print row_format.format(ds, s[ds]["AODstatus"], s[ds]["DAODstatus"], dsName)

def main():

    parser = argparse.ArgumentParser(description='Check status of datasets based on DSIDs and AMI tags. Written by C. Ohm & M. Tripiana')
    parser.add_argument('-p', '--projectTag', type=str, nargs='?', help='Project tag, defaults to "mc15_13TeV"', default='mc15_13TeV')
    parser.add_argument('-d', '--dsids', type=str, nargs='?', help='Text file(s) containing DSIDs', default='')
    parser.add_argument('-f', '--format', type=str, help='Format: DAOD_SUSY1 (default), DAOD_SUSY2, ...', default='DAOD_SUSY1')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Verbose mode, more detailed output and commands, etc')
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'), help='Save the dict holding the sample info to a pickle file')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), help='Open a pickle file containing a dict from a previous session')

    args = parser.parse_args()
    if args.verbose:
        print args

    if args.infile:
        with args.infile as handle:
            samples = pickle.load(handle)        
            printSamplesDict(samples)
            sys.exit()

    if not args.dsids:
        parser.print_help()
        sys.exit()

    f = open(args.dsids)

    dsids = [line.rstrip('\n') for line in f]

    if args.projectTag.startswith('data'):
        dsids = [ '00%s' % i for i in dsids if i and not i.startswith('#') ]

    samples = OrderedDict()

    for dsid in dsids:
        if not dsid or dsid.startswith('#'):
            continue
        #if args.verbose:
        print "Checking DSID %s..." % dsid
        samples[dsid] = getSamplesFromPattern(dsid, args.projectTag, defaultSimTags, defaultRecoTags, args.verbose, args.format, defaultDerivationTags)

    if args.verbose:
        print "These samples were found:"
        for dsid in samples:
            print "DSID: %s" % dsid
            for ds in samples[dsid]:
                print "  %s: %s" % (ds, samples[dsid][ds])

    # now check the status of each according to AMI:
    for dsid in samples:
        for ds in samples[dsid].keys():
            if ds == " N/A ":
                samples[dsid][ds+"status"] = " N/A "
            samples[dsid][ds+"status"] = getAmiStatus(samples[dsid][ds], args.verbose)

    print "Done, here is the status of your samples"
    if args.verbose:
        print "Here's the dict holding all the extracted info:"
        print samples

    if args.outfile:
        with args.outfile as handle:
            pickle.dump(samples, handle)
        print "Python dict with sample info saved to %s - you can read it in with the -i option" % args.outfile.name

    printSamplesDict(samples)

if __name__ == '__main__':
    main()
