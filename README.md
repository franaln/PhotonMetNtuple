PhotonMetNtuple
===============

## Current Base and SUSYTools release

    AnalysisBase,2.4.21 + SUSYTools-00-08-14


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


## Run

* To test locally:

    ```
    run.py --test /path/to/samples
    ```

* And to run in the grid:

    ```
    run.py -i input.txt [-d DIDS] --grid -v XX
    ```



