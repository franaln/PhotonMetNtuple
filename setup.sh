export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

##localSetupDQ2Client --skipConfirm
##localSetupFAX --skipConfirm
localSetupPandaClient ##currentJedi --noAthenaCheck

#rcSetup SUSY,2.3.23
#rcSetup SUSY,2.3.32
rcSetup SUSY,2.3.33a 

export PATH="$PWD/scripts:$PATH"
