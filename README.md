PhotonMetNtuple
===============

## Current Base and SUSYTools release

    AnalysisBase,2.4.24 (+ SUSYTools-00-08-33)


## Compile

    source setup.sh    

    rc checkout packages.txt
    rc checkout SUSYTools/doc/packages.txt

    rc find_packages
    rc clean
    rc compile


## ROOT Files needed

* PRW MC and ilumicalc files used from cvmfs using PathResolver, still configured via configfile


## Run

* To test locally:

    ```
    run.py --test /path/to/samples
    ```

* And to run in the grid:

    ```
    run.py -i input.txt [-d DIDS] --grid -v XX
    ```



