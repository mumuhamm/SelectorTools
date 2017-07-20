#!/bin/bash

lumi=35.87
lepid="Tight"
DATE=$(date +%d%b%Y)
#frfile=/eos/user/k/kelong/WZAnalysisData/FakeRates/fakeRate${DATE}-3LooseLeptons-${lepid}Leps.root
frfile=/eos/user/k/kelong/WZAnalysisData/FakeRates/fakeRate26Jun2017-3LooseLeptons-${lepid}Leps.root
input=WselectionLooseLeps

for jPt1 in `seq 35 5 50`; do
    for jPt2 in `seq 30 5 $jPt1`; do
        output=VBSselection_j1_${jPt1}_j2_${jPt2}
        histfile=/eos/user/k/kelong/WZAnalysisData/HistFiles/${output}-${DATE}-${lepid}.root

        cd $CMSSW_BASE/src/Analysis/WZAnalysis
        #./Utilities/scripts/makeFakeRates.py -s 3LooseLeptons -l $lumi -o $frfile
#        python ScaleFactors/setupScaleFactors.py -t $frfile 
#        echo ./Utilities/scripts/makeHistFile.py -l $lumi \
#            -s $input -o $histfile --output_selection $output \
#            -b "mjj,mjj_jesUp,mjj_jesDown,mjj_jerUp,mjj_jerDown,jetPt[0],jetPt[1],jetPt[2],jetEta[0],jetEta[1],jetEta[2],dEtajj,nJets"

        echo ./Utilities/scripts/prepareCombine.py \
            --input_file $histfile --output_file /eos/user/k/kelong/WZAnalysisData/CombineData/VBSselection_JetPtOptimization/`basename $histfile` \
            -l $lumi -s WZxsec2016/VBSselection \
            --folder_name j1Pt${jPt1}-j2Pt${jPt2}/$DATE 
    done
done

