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

    lines = open(args.input_file).read().split('\n')

    for line in lines:
        if not line or line.startswith('#'):
            continue

        if dids:
            did = int(line.split('.')[1])
            if did not in dids:
                continue

        samples.append(line)

    return samples



def run_job(sample, driver):

    if args.output is None:
        args.output = 'output'
    shutil.rmtree(args.output, True)

    short_name = '.'.join(sample.split('.')[0:3])

    # create a new sample handler to describe the data files we use
    logging.info("creating new sample handler")
    sh = ROOT.SH.SampleHandler()

    if driver == 'grid':
        ROOT.SH.scanDQ2(sh, sample)
    else:
        ROOT.SH.scanDir(sh, sample)

    # set the name of the tree in our files
    sh.setMetaString("nc_tree", "CollectionTree")

    logging.info("creating new job")
    job = ROOT.EL.Job()
    job.sampleHandler(sh)

    if args.nevents:
        logging.info("processing only %d events", args.nevents)
        job.options().setDouble(ROOT.EL.Job.optMaxEvents, args.nevents)
    
    # add our algorithm to the job
    logging.info("creating algorithms")
    
    if job_name == 'xAODAnalysis':

        if args.config_file is None:
            logging.error('you need to provide a configfile!')

        alg = ROOT.xAODAnalysis()

        is_data     = ('data15' in sample)
        is_susy     = ('_GGM_' in sample)
        is_susy_ewk = ('_GGM_mu_' in sample)
        is_atlfast  = (is_susy or 'MadGraphPythia8EvtGen_A14NNPDF23LO_ttgamma' in sample)
        
        logging.info('config: is data     = %s' % is_data)
        logging.info('config: is MC       = %s' % (not is_data))
        logging.info('config: is susy     = %s' % is_susy)
        logging.info('config: is susy ewk = %s' % is_susy_ewk)
        logging.info('config: is atlfast  = %s' % is_atlfast)
        logging.info('config: dosyst      = %s' % args.dosyst)


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
        
        #if job_name == 'xAODAnalysis':
        outname = 'user.' + os.environ['USER'] + '.' + short_name + '.mini_v' + args.version
        # elif job_name == 'xAODCountEwkProcesses':
        #     outname = 'user.' + os.environ['USER'] + '.' + short_name + '.ewk_v' + args.version

        driver.options().setString('nc_outputSampleName', outname)
        driver.options().setString(ROOT.EL.Job.optGridExcludedSite, 'ANALY_RHUL_SL6,ANALY_QMUL_SL6,ANALY_QMUL_HIMEM_SL6')
        driver.options().setString(ROOT.EL.Job.optGridNGBPerJob, 'MAX')
        driver.options().setString(ROOT.EL.Job.optGridMergeOutput, 'FALSE')

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
    parser.add_argument("--testdata", action='store_true')
    parser.add_argument("--testmc", action='store_true')
    parser.add_argument("--testsig", action='store_true')
    parser.add_argument('--test', dest='test')

    # run
    parser.add_argument('-i', dest='input_file')
    parser.add_argument('-s', dest='samples')
    parser.add_argument('-d', dest='dids', type=str)

    parser.add_argument('--job', default='xAODAnalysis')

    parser.add_argument('-c', dest='config_file', help='Config file')
    parser.add_argument('-v', dest='version')

    parser.add_argument("--dosyst", action='store_true', help="Create systematics blocks")

    parser.add_argument('--grid', action='store_true')    
    parser.add_argument('--dry', action='store_true')

    # others
    parser.add_argument('--download', action='store_true')


    
    global args
    args = parser.parse_args()

    if args.testdata or args.testmc or args.testsig or args.test is not None:
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


    logging.info("loading packages")
    ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

    global job_name
    job_name = args.job

    # if job_name == 'xAODAnalysis':

        # direc = os.environ["ROOTCOREBIN"] + "/data/PhotonMetNtuple/"
        # os.system('mkdir -p %s' % direc)

        # # GRL
        # #grl = '/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml'
        # grl = '/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml' # new 20.7

        # shutil.copy2(grl, direc + "grl.xml")

        # # ST config
        # conf = 'ST_PhotonMet.conf'
        # shutil.copy2('data/'+conf, direc + conf)

        # # PU files
        # pu_files = [
        #     'merged_prw_mc15c.root',
        #     'ilumicalc_histograms_None_276262-284154.root',
        #     # 'signal_prw.root',
        #     # 'signal_ewk_prw.root',
        #     # 'merged_prw.root',
        #     # 'ttgam_prw.root',
        #     ]

        # for f in pu_files:
        #     shutil.copy2('/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/'+f, direc+f)


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


    # Test code
    elif args.testmc:
        args.output = 'output_mc'
        run_job('/ar/pcunlp001/raid/falonso/SUSY1/mc_gamjet', 'local')
        
    elif args.testdata:
        args.output = 'output_data'
        run_job('/ar/pcunlp001/raid/falonso/SUSY1/data', 'local')

    elif args.testsig:
        args.output = 'output_signal'
        run_job('/ar/pcunlp001/raid/falonso/SUSY1/mc15_13TeV.373061.HerwigppEvtGen_UEEE5CTEQ6L1_GGM_M3_mu_1100_450.merge.AOD.e4349_a766_a777_r6282', 'local')

    elif args.test is not None:
        if args.output is None:
            args.output = 'output_test' 
        run_job(args.test, 'local')


if __name__ == "__main__":
    main()
