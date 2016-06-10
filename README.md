PhotonMetNtuple
===============

## Recommended ST tag for 20.7 (probably this pkg is not compatible with 20.1 anymore, use latest tag)

    Base,2.4.10 + SUSYTools-00-07-65

## To compile:

    source setup.sh    

    rc checkout_pkg $SVNOFF/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/SUSYTools-00-07-65
    rc checkout SUSYTools/doc/packages.txt

    rc find_packages
    rc clean
    rc compile


## Run

    run.py -i input.txt --grid -v XX

