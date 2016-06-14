export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

localSetupPandaClient

rcSetup Base,2.4.11  # and checkout SUSYTools-00-07-69

export PATH="$PWD/scripts:$PATH"
