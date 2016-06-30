export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

localSetupPandaClient

rcSetup Base,2.4.13 # + SUSYTools-00-07-82

export ROOTCORE_NCPUS="4"

export PATH="$PWD/scripts:$PATH"
