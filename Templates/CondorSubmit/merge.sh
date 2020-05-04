#!/bin/bash
CMSSW_RELEASE_BASE="${CMSSW_RELEASE_BASE}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd $$CMSSW_RELEASE_BASE
eval `scramv1 runtime -sh`
popd
export LD_LIBRARY_PATH=$$PWD/lib:$$LD_LIBRARY_PATH

outfile=$$1
if [[ $$1 == /eos/user* ]]; then
    outfile="root://eosuser.cern.ch/$$1"
fi

hadd -f $$outfile *.root
