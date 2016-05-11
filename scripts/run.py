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


def run_job(sample, driver):

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
        alg = ROOT.xAODAnalysis()

        is_data    = ('data' in sample)
        is_atlfast = ('_GGM_' in sample)

        print 'data:   ', is_data
        print 'atlfast:', is_atlfast
    
        alg.isData = is_data
        alg.isAtlfast = is_atlfast
        alg.doSyst = args.dosyst

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
        driver.options().setString(ROOT.EL.Job.optGridExcludedSite, 'CA ANALY_TRIUMF,US ANALY_HU_ATLAS_Tier2,CA ANALY_SFU,FR ANALY_GRIF-IRFU')
        driver.options().setString(ROOT.EL.Job.optGridNGBPerJob, 'MAX')
        driver.options().setString(ROOT.EL.Job.optGridMergeOutput, 'FALSE')

        logging.info('submit job: ' + outname)
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
    parser.add_argument("--output", help="dir to store the output", default="output")
    #parser.add_argument("--inputDS", help="input DS from DQ2")
    #parser.add_argument("--driver", help="select where to run", choices=("direct", "prooflite", "grid"), default="direct")
    parser.add_argument("--nevents", type=int, help="number of events to process for all the datasets")
    # parser.add_argument("--skip-events", type=int, help="skip the first n events")

    # parser.add_argument("--isdata", action='store_true')
    # parser.add_argument("--isatlfast", action='store_true')
    parser.add_argument("--dosyst", action='store_true', help="Create Trees with systemtic variations")
    
    parser.add_argument("--testdata", action='store_true')
    parser.add_argument("--testmc", action='store_true')
    parser.add_argument("--testsig", action='store_true')
    parser.add_argument('--test', dest='test')

    parser.add_argument('--grid', action='store_true')    
    parser.add_argument('-v', dest='version')
    parser.add_argument('-s', dest='samples')
    parser.add_argument('-d', dest='dry')

    parser.add_argument('--download', action='store_true')

    # job
    parser.add_argument('--job', default='xAODAnalysis')


    
    global args
    args = parser.parse_args()

    if args.testdata or args.testmc or args.testsig or args.test is not None:
        args.version = 0 

    if args.version is None:
        parser.print_usage()
        return 1


    if args.download:

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

    if job_name == 'xAODAnalysis':
    
        grl = '/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml'

        pu_files = [
            'ilumicalc_histograms_None_276262-282712.root',
            'signal_prw.root',
            'signal_ewk_prw.root',
            'merged_prw.root',
            'ttgam_prw.root',
            'SUSYTools_Default.conf',
            ]

        direc = os.environ["ROOTCOREBIN"] + "/data/PhotonMetNtuple/"
        newname = direc + "grl.xml"
        shutil.copy2(grl, newname)

        for f in pu_files:
            shutil.copy2('data/'+f, direc+f)


    #Set up the job for xAOD access:
    ROOT.xAOD.Init().ignore()

    if args.grid:

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
        run_job('/raid/falonso/SUSY1/mc_gamjet', 'local')
        
    elif args.testdata:
        args.output = 'output_data'
        run_job('/raid/falonso/SUSY1/data', 'local')

    elif args.testsig:
        args.output = 'output_signal'
        run_job('/raid/falonso/signal', 'local')

    elif args.test is not None:
        args.output = 'output_test'
        run_job(args.test, 'local')


if __name__ == "__main__":
    main()
