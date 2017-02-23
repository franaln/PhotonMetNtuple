#! /usr/bin/env python

import os
import ROOT
import logging
import shutil
import argparse
import subprocess

ROOT.PyConfig.IgnoreCommandLineOptions = True

logging.basicConfig(level=logging.INFO)

excluded_sites = 'ANALY_RHUL_SL6' #ANALY_CONNECT_SHORT,ANALY_INFN-ROMA1,ANALY_wuppertalprod' #ANALY_RHUL_SL6,ANALY_QMUL_SL6,ANALY_QMUL_HIMEM_SL6,ANALY_AGLT2_SL6,ANALY_RAL_SL6,'

import atexit
@atexit.register
def quite_exit():
    ROOT.gSystem.Exit(0)

def get_grid_name(sample, version):

    splitted_sample = sample.split('.')

    short_name = '.'.join(splitted_sample[:3])

    tags = splitted_sample[-1]
    ptag = tags.split('_')[-1] if tags.split('_')[-1].startswith('p') else ''

    if ptag:
        outname = 'user.' + os.environ['USER'] + '.' + short_name + '.mini.' + ptag + '.v' + args.version
    else:
        outname = 'user.' + os.environ['USER'] + '.' + short_name + '.mini.v' + args.version

    return outname

def check_grid_exists(sample, version):

    grid_name = get_grid_name(sample, version) + '_output.root'

    cmd = 'rucio ls --short --filter type=CONTAINER %s | sort -r ' % grid_name
    cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    query_result = cmd_output.stdout

    lines = [line.rstrip() for line in query_result.readlines()]
    print lines
    return True if lines else False


def get_samples_from_file(file_, dids_str=None):

    logging.info('get samples from %s with dids = %s' % (file_, dids_str))

    include_dids = []
    exclude_dids = []
    if dids_str is not None:
        try:
            if '-' in dids_str:
                first, last = dids_str.split('-')
                if first:
                    include_dids = range(int(first), int(last)+1)
                else:
                    exclude_dids = [int(last), ]

            elif ',' in dids_str:
                include_dids = [ int(did) for did in dids_str.split(',') ]
            else:
                include_dids = [int(dids_str),]
                
        except:
            logging.error('bad dids syntax. ignoring...')

    samples = []

    for line in open(file_).read().split('\n'):
        if not line or line.startswith('#'):
            continue

        if include_dids or exclude_dids:
            did = int(line.split('.')[1])

            if include_dids and did not in include_dids:
                continue
            if exclude_dids and did in exclude_dids:
                continue

        samples.append(line.strip())

    if not samples:
        logging.error('not samples to run')

    # for s in samples:
    #     print check_grid_exists(s, args.version)

    return samples


def run_job(sample, driver):

    logging.info('running job (%s) for sample: %s' % (driver, sample))

    if args.output is None:
        args.output = 'output'
    shutil.rmtree(args.output, True)

    # create a new sample handler to describe the data files we use
    sh = ROOT.SH.SampleHandler()
    if driver == 'grid':
        ROOT.SH.scanRucio(sh, sample)
    else:
        ROOT.SH.scanDir(sh, sample)

    # set the name of the tree in our files
    sh.setMetaString("nc_tree", "CollectionTree")

    job = ROOT.EL.Job()
    job.sampleHandler(sh)

    job.options().setDouble(ROOT.EL.Job.optXAODSummaryReport, False)

    if args.nevents:
        logging.info("processing only %d events", args.nevents)
        job.options().setDouble(ROOT.EL.Job.optMaxEvents, args.nevents)
    
    # add our algorithm to the job
    if alg_name in ['xAODAnalysis', 'xAODBaselineAnalysis', 'xAODJfakeSample']:

        alg = getattr(ROOT, alg_name)()
        print alg
        alg.mem = [] # use this to prevent ownwership problems
        
        is_data     = ('data15' in sample or 'data16' in sample)
        is_susy     = ('_GGM' in sample)
        is_susy_ewk = ('_GGM_mu' in sample)
        is_atlfast  = (is_susy or 'MadGraphPythia8EvtGen_A14NNPDF23LO_ttgamma' in sample)
        

        logging.info('--')
        logging.info('-- Configuration to use. Please check if it is ok!')
        logging.info('-- Configfile  = %s' % args.config_file)
        logging.info('-- Sample = %s' % sample)
        logging.info('-- is data=%s, is MC=%s, is atlfast=%s, is susy=%s, is susy ewk=%s' % (is_data, (not is_data), is_atlfast, is_susy, is_susy_ewk))
        logging.info('-- do systematics = %s' % args.dosyst)
        logging.info('--')

        alg.config_file = args.config_file
        alg.is_data = is_data
        alg.is_susy = is_susy
        alg.is_susy_ewk = is_susy_ewk
        alg.is_atlfast = is_atlfast
        alg.do_syst = args.dosyst

    elif alg_name == 'xAODCountEwkProcesses':
        alg = ROOT.xAODCountEwkProcesses()

    elif alg_name == 'xAODTruthAnalysis':
        alg = ROOT.xAODTruthAnalysis()
        if args.dopdfrw:
            alg.do_pdfrw = True


    logging.info("adding algorithms")
    job.algsAdd(alg)

    # make the driver we want to use:
    # this one works by running the algorithm directly
    if driver == 'grid':

        logging.info('running on Grid') 

        driver = ROOT.EL.PrunDriver() 
        
        # output name
        splitted_sample = sample.split('.')

        short_name = '.'.join(splitted_sample[0:3])

        if alg_name == 'xAODAnalysis':
            outname = get_grid_name(sample, args.version)

        elif alg_name == 'xAODBaselineAnalysis':
            outname = 'user.' + os.environ['USER'] + '.' + short_name + '.base.v' + args.version

        elif alg_name == 'xAODCountEwkProcesses':
            outname = 'user.' + os.environ['USER'] + '.' + short_name + '.ewk.v' + args.version

        elif alg_name == 'xAODTruthAnalysis':
            outname = 'user.' + os.environ['USER'] + '.' + short_name + '.truth.v' + args.version

        # driver options
        driver.options().setString('nc_outputSampleName', outname)
        if excluded_sites:
            driver.options().setString(ROOT.EL.Job.optGridExcludedSite, excluded_sites)
        driver.options().setString(ROOT.EL.Job.optGridNGBPerJob, 'MAX')
        driver.options().setString(ROOT.EL.Job.optGridMergeOutput, 'FALSE')
        driver.options().setDouble(ROOT.EL.Job.optRemoveSubmitDir, 1)

        # submit 
        logging.info('submit job: ' + outname)
        if not args.dry:
            driver.submitOnly(job, args.output)

    # elif driver == 'batch':

    #     logging.info('running on prooflite')
    #     driver = ROOT.EL.ProofDriver()
    #     logging.info('submit job')
    #     driver.submit(job, args.output)
        
    elif driver == 'local':

        logging.info('running on direct')
        driver = ROOT.EL.DirectDriver()
        logging.info('submit job')

        try:
            driver.submit(job, args.output)
        except KeyboardInterrupt:
            raise


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="dir to store the output")
    parser.add_argument("--nevents", type=int, help="number of events to process for all the datasets")

    parser.add_argument('--alg', default='xAODAnalysis')
    parser.add_argument("--dosyst", action='store_true', help="Create systematics blocks")

    # test
    parser.add_argument('--test', dest='test')

    # run (in the grid)
    parser.add_argument('--grid', action='store_true')    
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--dry', action='store_true')

    parser.add_argument('-i', dest='input_file')
    parser.add_argument('-d', dest='dids', type=str)

    parser.add_argument('-c', dest='config_file', default='PhotonMetNtuple_20.7_std.conf', help='Config file')
    parser.add_argument('-v', dest='version')

    parser.add_argument('--dopdfrw', action='store_true')
    # parser.add_argument('--dotar', dest='do_tarball', action='store_true')
    # parser.add_argument('--usetar', dest='use_tarball', action='store_true')
    
    global args
    args = parser.parse_args()

    if args.show:

        if args.input_file is not None:
            torun = get_samples_from_file(args.input_file, args.dids)

        for sample in torun:
            print sample
        
        return 0


    if args.test is not None:
        args.version = 0 

    if args.version is None:
        parser.print_usage()
        return 1

    if args.download:

        if args.input_file is not None:
            torun = get_samples_from_file(args.input_file, args.dids)

        print 'rucio download \\'
        for i, sample in enumerate(torun):
            
            outname = get_grid_name(sample, args.version) + '_output.root'

            if i == len(torun) - 1:
                print '    %s' % outname
            else:
                print '    %s \\' % outname

        return 0


    global alg_name
    alg_name = args.alg

    if alg_name == 'xAODAnalysis' and args.config_file is None:
        logging.error('you need to provide a configfile!')

    ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
    ROOT.xAOD.Init().ignore()

    # Run on the GRID
    if args.grid:
        if args.input_file is not None:
            torun = get_samples_from_file(args.input_file, args.dids)

        for sample in torun:
            run_job(sample, 'grid')


    # Test code locally
    elif args.test is not None:
        if args.output is None:
            args.output = 'output_test' 
        run_job(args.test, 'local')


if __name__ == "__main__":
    main()
