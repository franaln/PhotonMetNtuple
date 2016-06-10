PhotonMetNtuple
===============

## Recommended ST tag for 20.7 

    Base,2.4.10 + SUSYTools-00-07-65

* probably this pkg is not compatible with 20.1 anymore, use latest tag


## To compile:

    source setup.sh    

    rc checkout_pkg $SVNOFF/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/SUSYTools-00-07-65
    rc checkout SUSYTools/doc/packages.txt

    rc find_packages
    rc clean
    rc compile


## Run

* First you need copy the needed files (PRW files) to PhotonMetNtuple/data

* To test:

    ```
    run.py --test /path/to/samples -c PhotonMetNtuple_20.7_std.conf
    ```

* And to run in the grid:

    ```
    run.py -i input.txt --grid -v XX -c PhotonMetNtuple_20.7_std.conf
    ```

