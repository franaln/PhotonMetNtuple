PhotonMetNtuple
===============

## Current Base and SUSYTools release

    AnalysisBase,2.4.29 (+ SUSYTools-00-08-58)


## Compile

    source setup.sh    

    rc checkout packages.txt
    rc checkout SUSYTools/doc/packages.txt

    rc find_packages
    rc clean
    rc compile

## Make a new tag

    git tag -a vXX -m "tagging version XX"


## ROOT Files needed

* PRW MC and ilumicalc files used from cvmfs using PathResolver, still configured via configfile

* Yo may need to copy the data files from /afs/cern.ch/user/f/falonso/public/SUSY/PMN_data to PhotonMetNtuple/data/


## Run

* To test locally:

    ```
    run.py --test /path/to/samples
    ```

* And to run in the grid:

    ```
    run.py -i input.txt [-d DIDS] --grid -v XX
    ```
    
* To run the truth analysis (or any other analysis) use the --alg option:

    ```
    run.py --alg xAODTruthAnalysis ....
    ```
    