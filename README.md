PhotonMetNtuple
===============

## Current Base and SUSYTools release (see versions.md for more)

    AnalysisBase,2.4.19 + SUSYTools-00-08-04 + packages.txt


## Compile

    source setup.sh    

    rc checkout_pkg $SVNOFF/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/SUSYTools-00-XX-XX

    rc checkout SUSYTools/doc/packages.txt

    rc find_packages
    rc clean
    rc compile

## Files needed

* PRW MC files: /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/

* ilumicalc files: from Luminosity Calculator using corresponding GRL



## Run

* To test:

    ```
    run.py --test /path/to/samples -c PhotonMetNtuple_20.7_std.conf
    ```

* And to run in the grid:

    ```
    run.py -i input.txt -d [DIDS] --grid -v XX -c PhotonMetNtuple_20.7_std.conf
    ```

## Samples


