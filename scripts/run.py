#! /usr/bin/env python

import os
import ROOT
import logging
import shutil
import argparse
import samples

logging.basicConfig(level=logging.INFO)

import atexit
@atexit.register
def quite_exit():
    ROOT.gSystem.Exit(0)


def get_samples_from_file(file_, dids_str=None):

    logging.info('get samples from %s with dids = %s' % (file_, dids_str))

    dids = []
    if dids_str is not None:

        try:
            if '-' in dids_str:
                first, last = dids_str.split('-')
                dids = range(int(first), int(last)+1)
            elif ',' in dids_str:
                dids = [ int(did) for did in dids_str.split(',') ]
            else:
                dids = [int(dids_str),]
        except:
            logging.error('bad dids syntax. ignoring...')


    samples = []

    lines = open(file_).read().split('\n')

    for line in lines:
        if not line or line.startswith('#'):
            continue

        if dids:
            did = int(line.split('.')[1])
            if did not in dids:
                continue

        samples.append(line)

    if not samples:
        logging.error('not samples to run')

    return samples



def run_job(sample, driver):

    logging.info('running job (%s) for sample: %s' % (driver, sample))

    if args.output is None:
        args.output = 'output'
    shutil.rmtree(args.output, True)


    # create a new sample handler to describe the data files we use
    sh = ROOT.SH.SampleHandler()

    if driver == 'grid':
        ROOT.SH.scanDQ2(sh, sample)
    else:
        ROOT.SH.scanDir(sh, sample)

    # set the name of the tree in our files
    sh.setMetaString("nc_tree", "CollectionTree")

    job = ROOT.EL.Job()
    job.sampleHandler(sh)

    if args.nevents:
        logging.info("processing only %d events", args.nevents)
        job.options().setDouble(ROOT.EL.Job.optMaxEvents, args.nevents)
    
    # add our algorithm to the job
    if job_name == 'xAODAnalysis':

        alg = ROOT.xAODAnalysis()

        is_data     = ('data15' in sample)
        is_susy     = ('_GGM_' in sample)
        is_susy_ewk = ('_GGM_mu_' in sample)
        is_atlfast  = (is_susy or 'MadGraphPythia8EvtGen_A14NNPDF23LO_ttgamma' in sample)
        
        logging.info('-- Configuration to use. Please check if it is ok! --')
        logging.info('configfile  = %s' % args.config_file)
        logging.info('is data=%s, is MC=%s, is atlfast=%s, is susy=%s, is susy ewk=%s' % (is_data, (not is_data), is_susy, is_susy_ewk, is_atlfast))
        logging.info('do systematics = %s' % args.dosyst)
        logging.info('----')

        alg.config_file = args.config_file
        alg.is_data = is_data
        alg.is_susy = is_susy
        alg.is_susy_ewk = is_susy_ewk
        alg.is_atlfast = is_atlfast
        alg.do_syst = args.dosyst

    elif job_name == 'xAODCountEwkProcesses':
        alg = ROOT.xAODCountEwkProcesses()


    logging.info("adding algorithms")
    job.algsAdd(alg)

    # make the driver we want to use:
    # this one works by running the algorithm directly
    logging.info("creating driver")

    if driver == 'grid':

        logging.info('running on Grid') 

        driver = ROOT.EL.PrunDriver() 
        

        # output name
        splitted_sample = sample.split('.')

        short_name = '.'.join(splitted_sample[0:3])

        tags = splitted_sample[-1]
        ptag = tags.split('_')[-1] if tags.split('_')[-1].startswith('p') else ''

        if job_name == 'xAODAnalysis':
            if ptag:
                outname = 'user.' + os.environ['USER'] + '.' + short_name + '.mini.' + ptag + '.v' + args.version
            else:
                outname = 'user.' + os.environ['USER'] + '.' + short_name + '.mini.v' + args.version

        elif job_name == 'xAODCountEwkProcesses':
            outname = 'user.' + os.environ['USER'] + '.' + short_name + '.ewk_v' + args.version

        # driver options
        driver.options().setString('nc_outputSampleName', outname)
        driver.options().setString(ROOT.EL.Job.optGridExcludedSite, 'ANALY_RHUL_SL6,ANALY_QMUL_SL6,ANALY_QMUL_HIMEM_SL6')
        driver.options().setString(ROOT.EL.Job.optGridNGBPerJob, 'MAX')
        driver.options().setString(ROOT.EL.Job.optGridMergeOutput, 'FALSE')

        # submit 
        logging.info('submit job: ' + outname)
        if not args.dry:
            driver.submitOnly(job, args.output)

    elif driver == 'batch':

        logging.info('running on prooflite')
        driver = ROOT.EL.ProofDriver()
        logging.info('submit job')
        driver.submit(job, args.output)
        
    elif driver == 'local':

        logging.info('running on direct')
        driver = ROOT.EL.DirectDriver()
        logging.info('submit job')
        driver.submit(job, args.output)
        

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="dir to store the output")
    parser.add_argument("--nevents", type=int, help="number of events to process for all the datasets")

    # test
    parser.add_argument('--test', dest='test')

    # run
    parser.add_argument('-i', dest='input_file')
    parser.add_argument('-s', dest='samples')
    parser.add_argument('-d', dest='dids', type=str)

    parser.add_argument('--job', default='xAODAnalysis')

    parser.add_argument('-c', dest='config_file', default='PhotonMetNtuple_20.7_std.conf', help='Config file')
    parser.add_argument('-v', dest='version')

    parser.add_argument("--dosyst", action='store_true', help="Create systematics blocks")

    parser.add_argument('--grid', action='store_true')    
    parser.add_argument('--dry', action='store_true')

    # others
    parser.add_argument('--download', action='store_true')
    
    global args
    args = parser.parse_args()

    if args.test is not None:
        args.version = 0 

    if args.version is None:
        parser.print_usage()
        return 1


    if args.download:

        if args.input_file is not None:
            torun = get_samples_from_file(args.input_file, args.dids)

        elif args.samples is not None:

            try:
                torun = []
                for s in args.samples.split(','):
                    torun.extend(getattr(samples, s))
            except:
                print 'no samples to run'
                torun = []

        for sample in torun:
            s = '.'.join(sample.split('.')[:3])
            outname = 'user.falonso.%s.mini_v%s_output.root/' % (s, args.version)
            print 'rucio download %s' % outname

        return 0


    ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

    global job_name
    job_name = args.job

    if job_name == 'xAODAnalysis' and args.config_file is None:
        logging.error('you need to provide a configfile!')


    #Set up the job for xAOD access:
    ROOT.xAOD.Init().ignore()

    if args.grid:
        
        if args.input_file is not None:
            torun = get_samples_from_file(args.input_file, args.dids)

        elif args.samples is not None:

            try:
                torun = []
                for s in args.samples.split(','):
                    torun.extend(getattr(samples, s))
            except:
                raise

        for sample in torun:
            run_job(sample, 'grid')


    # Test code locally
    elif args.test is not None:
        if args.output is None:
            args.output = 'output_test' 
        run_job(args.test, 'local')


if __name__ == "__main__":
    main()
