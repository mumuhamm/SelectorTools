#!/bin/bash

dataset_manager=$(./Utilities/scripts/getConfigValue.py dataset_manager_path)/AnalysisDatasetManager
pushd $CMSSW_BASE/src/Analysis/SelectorTools/Cuts

echo "INFO: Linking alias files"
for folder in $(ls -d */); do
    analysis_name=${folder/\//}
    alias_file=${dataset_manager}/Aliases/${analysis_name}.json 
    if [ -f $alias_file ]; then
        pushd $analysis_name
        ln -s $alias_file aliases.json
        echo "INFO: Linked analysis $analysis_name to alias file"
        popd
    fi
done

popd

analysis=$(./Utilities/scripts/getConfigValue.py analysis)
if [[ $analysis == WZ ]]; then
    echo "INFO: Downloading scale factor files"
    pushd ScaleFactors
    bash setup.sh
    popd
fi

if [ -f PileupWeights/calculatePileupCorrections.py ]; then
    echo "INFO: Producing pileup weight distribution files"
    pushd PileupWeights
    bash getDataDistribution.sh
    python calculatePileupCorrections.py
    popd
else
    echo "WARNING! PileupWeights repository not found. You can clone with"
    echo ' --recursive, or run "git submodule update --init" now'
    echo 'NOTE: Not needed for gen-only studies"
fi
