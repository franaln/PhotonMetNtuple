PhotonMetNtuple
===============

## Current Base and SUSYTools release

    AnalysisBase,2.4.22 (+ SUSYTools-00-08-24)


## Compile

    source setup.sh    

    rc checkout packages.txt
    rc checkout SUSYTools/doc/packages.txt

    rc find_packages
    rc clean
    rc compile


## Files needed

* PRW MC files: /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/

* ilumicalc files: from Luminosity Calculator using corresponding GRL

* OR copy all needed files from /afs/cern.ch/user/f/falonso/public/SUSY/PMN_data to PhotonMetNtuple/data/


## Run

* To test locally:

    ```
    run.py --test /path/to/samples
    ```

* And to run in the grid:

    ```
    run.py -i input.txt [-d DIDS] --grid -v XX
    ```



