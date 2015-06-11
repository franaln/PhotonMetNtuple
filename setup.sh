export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
localSetupDQ2Client --skipConfirm
localSetupFAX --skipConfirm
localSetupPandaClient currentJedi --noAthenaCheck
rcSetup Base,2.3.12
# rc find_packages
# rc compile