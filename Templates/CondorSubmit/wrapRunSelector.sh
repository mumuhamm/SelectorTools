#!/bin/bash

# Following implementation by N. Smith, Fermilab
# https://gitlab.cern.ch/ncsmith/monoZ/tree/master/selector 

# TODO: should be filled in from a template
tar xvzf ${tarball}
CMSSW_RELEASE_BASE="${CMSSW_RELEASE_BASE}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd $$CMSSW_RELEASE_BASE
eval `scramv1 runtime -sh`
popd
export LD_LIBRARY_PATH=$$PWD/lib:$$LD_LIBRARY_PATH
./Utilities/scripts/makeHistFile.py $$@ || exit $$?

