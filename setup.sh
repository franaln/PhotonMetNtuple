export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

localSetupPandaClient

# rcSetup SUSY,2.3.45 # (SUSYTools 07-41)
# rcSetup SUSY,2.3.48 # (SUSYTools 07-45)
rcSetup Base,2.3.50  # SUSYTools-00-07-52)
export PATH="$PWD/scripts:$PATH"
