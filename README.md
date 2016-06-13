PhotonMetNtuple
===============

## Release

    Base,2.4.10 + SUSYTools-00-07-65 + METUtilities-00-02-29


## Compile

    source setup.sh    

    rc checkout_pkg $SVNOFF/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/SUSYTools-00-07-68
    rc checkout_pkg atlasoff/Reconstruction/MET/METUtilities/tags/METUtilities-00-02-29

    rc find_packages
    rc clean
    rc compile


## Run

* First you need copy the needed files (PRW files) to PhotonMetNtuple/data

    ```
    GRLs here: /afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/
    PRW here: 
    ```

* To test:

    ```
    run.py --test /path/to/samples -c PhotonMetNtuple_20.7_std.conf
    ```

* And to run in the grid:

    ```
    run.py -i input.txt --grid -v XX -c PhotonMetNtuple_20.7_std.conf
    ```

## Samples


## Data 

    data15: (20.7 repro): data15_13TeV.periodAllYear_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml, corresponding to ...
    data16: data16_13TeV.periodAllYear_DetStatus-v78-pro20-04_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml, corresponding to 2613.83 pb-1