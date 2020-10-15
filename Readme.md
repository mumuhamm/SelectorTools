Analysis code for WZ/ZZ analyses. Some scripts using selections to skim Ntuples and C++ code to make histogram files for WZ.

+ [Setup](#setup)
+ [Overview](#overview)
+ [Specifics](#overview)
    + [NanoAOD](#nanoaod)
        + [Skimming NanoAOD](#skimming-nanoaod)
     + [UWVV nutples](#producing-uwvv-ntuples)
        + [Skimming UWVV ntuples](#skimming-uwvv-ntuples)
    + [Producing processed histogram file](#running-analysis-code)
        + [Z selector example](#the-z-selector-example)
        + [Writing your own selector](#implementing-your-own-selector)
        + [WZ selector](#the-wz-selector)
        + [WZ Background estimation](#nonprompt-background-estimate-for-wz)
    + [Statistical analysis](#running-statistical-analysis)
    + [Producing plots](#plotting)
    

# Setup
-----------
CMSSW version: CMSSW_11_0_0
```bash
cmssw_version="11_0_0" # or anything that isn't too old
username="kdlong" # or your username
```

To checkout and compile:

```console
cmsrel CMSSW_version
cd CMSSW_${cmssw_version}/src
mkdir Analysis
cd Analysis
git clone git@github.com:<username>/SelectorTools.git
scram b -j 8
```

You will also want to install a separate package that contains information on the datasets used and histograms and binning. It is not required that you use this file, but it will be convenient for managing larger datasets (e.g., for getting the correct cross sections). It's recommended that you fork this into your github because some files/settings will be user specific.

install_path = ~username/work/ (or whatever you prefer)

```console
cd install_path
git clone git@github.com:<username>/AnalysisDatasetManager.git
```

You should create a new configuration file following the example [here](../master/Templates/config.kdlong) to have settings specific to you. Use the name config.<username>.
	
NOTE: Everytime you run the code, you need to be in this CMS environment, that is, run ```cmsenv``` from the CMSSW_X_Y_Z/src directory where the package is installed.

# Overview
-----------
This repository includes scripts to analyze ROOT ntuples using implementations of the [TSelector](https://root.cern.ch/doc/master/classTSelector.html) class. It was originally developed for [SMP-18-001](http://cms.cern.ch/iCMS/analysisadmin/cadilines?line=SMP-18-001), which is still reflected in the WZSelector classes. Because the structure is quite general, it has been extended for other CMS results, as well as in general analysis based on NanoAOD, at either the RECO or GEN level.

Most of the current development is devoted to analyses based on NaanoAOD, a common ntuple format in CMS. More detail is given in the [NanoAOD](#nanoaod) section.

A general analysis will proceed in several steps.

1. Produce ntuples. This step is independent of this package. If using NanoAOD, it is generally not necessary to produce your own samples. If using NanoGen, some useful example scripts are [here](https://github.com/kdlong/WMassNanoGen). For UWVV, see [UWVV](https://github.com/uwcms/UWVV). 
2. Skim ntuples or NanoAOD to create smaller files that can be copied to a "local" storage disk, such as /eos at CERN, /data at uwlogin or /nfs_scratch at Wisconsin. For NanoAOD, there are dedicated modules in the [NanoAODTools package](https://github.com/cms-nanoAOD/nanoAOD-tools).
3. Run analysis code to produce histograms with final selections.
4. Plotting code combines histograms by scaling with cross section and luminosity information. Colors and CMS style are implemented. Mostly outside the scope of this package.
5. Perform fit (using HiggsCombine package)

# Specifics
-----------
Each step deserves some degree of explanation. They are also all driven by independent scripts, and can be run separately. The ntuple skimming is generally run independently. It is not required, but is advantageous to reduce files sizes for convenience in storage and processing. It is absolutely necessary that the tightest condition you use in your skim be looser the selection you implement in a later selector, however. For a fully MC driven analysis, one can often implement the full selection at the skim step and produce plots from here. Note, however, that the statistical tools are designed to run over the output of the histogram files produced by the selector.

## NanoAOD

https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD

### Skimming NanoAOD

https://github.com/cms-nanoAOD/nanoAOD-tools

Dilepton and 3 lepton specific skimming tools:

https://github.com/kdlong/NanoVVSkims

## Running analysis code

The anlaysis code is driven by the script [makeHistFile.py](Utilities/scripts/makeHistFile.py). This script takes care of configuring the input and options and calling the C++ selector for the appropriate configuration with the appropriate arguments. Run ```./Utilities/scripts/makeHistfile.py --help``` to see the list of options.


### Gen analysis on NanoAOD

The Gen analysis on NanoAOD is based on the class [NanoGenSelectorBase.cc](src/NanoGenSelectorBase.cc). This class reads the variables from a NanoAOD and configures the job. Specific anlayses for W and Z selections are implemented in [ZGenSelector.cc](src/ZGenSelector.cc) and [WGenSelector.cc](src/WGenSelector.cc). You will need the [AnalysisDatasetManager repository](https://github.com/kdlong/AnalysisDatasetManager) to define the histogram binning. If you are only running the W and Z analyses, you can clone it using the sparse checkout step

```bash
bash <(curl -s https://raw.githubusercontent.com/kdlong/AnalysisDatasetManager/master/SparseCheckouts/sparseCheckout.sh) VGenSparse.config
```

You should then put the path to this repo in your Templates/user.config file inside the SelectorTools repository.

Test commands to run these analyses are

```./Utilities/scripts/makeHistFile.py -f /store/mc/RunIISummer16NanoAODv7/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/10000/AC80C98C-F0A6-7443-8DD3-68F5D09F0CAE.root -a ZGen -s None --input_tier NanoAOD -o testZ.root```

```./Utilities/scripts/makeHistFile.py -f /store/mc/RunIISummer16NanoAODv7/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/260000/AFEFB52A-AB45-334A-AFBD-C2FD0D28F3D6.root -a WGen -s None --input_tier NanoAOD -o testW.root```

Try these commands out and take a look at the output files (testZ.root and testW.root). Many histograms should be produced inside the main folder, which is name "Unknown," because you haven't specified a name. You can browse using the TBrowser, for example. You can specify a name for your folder using name@path, for example, to name the folder ZMiNNLO:

```./Utilities/scripts/makeHistFile.py -f ZMiNNLO@/store/mc/RunIISummer16NanoAODv7/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/10000/AC80C98C-F0A6-7443-8DD3-68F5D09F0CAE.root -a ZGen -s None --input_tier NanoAOD -o testZ.root```

Using this you can also process multiple files from different datasets:

```./Utilities/scripts/makeHistFile.py -f ZMGNLO@/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/120000/0643F6C0-C42B-564C-A9BD-161948BBFB81.root,ZMiNNLO@/store/mc/RunIISummer16NanoAODv7/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/10000/AC80C98C-F0A6-7443-8DD3-68F5D09F0CAE.root -o testZ.root```

Now you should have two folders in your file, one named ZMGNLO and one named ZMiNNLO. These are meant to be used for comparisons of the same distribution from a different dataset. You can play around in ROOT to overlay plots of the same quanity (e.g., ptZ_mm).

### The Z selector example

See [this commit](https://github.com/kdlong/VVAnalysis/commit/18a1d903e149653fff3985b43f1acc834632a7ac) for an example of how to add histograms to a selector using the configuration setup. [These corresponding changes](https://github.com/kdlong/AnalysisDatasetManager/commit/39909f6e76046b6ab39b293fa3ea209d1cd202a8) to the [AnalysisDatasetManager repository](https://github.com/kdlong/AnalysisDatasetManager) are necessary.

The main program that runs the analysis is the file ```./Utilities/scripts/makeHistFile.py``` . This file is a wrapper for the running the selector defined in your src directory over the files specified in the AnalysisDatasetManager. This code works with a conviention used in most of this code suite, that is a selection, analysis, and input tier are needed for most all functions to run. for completeness, each option for the code will be explained and the user can know how to use the software:

* **-s**: Selection. User defined selection to be implimented in Analyzer. Need to put in flags in cc files for make any difference (ie PlotObjects area)
* **-v**: Version. Pretty self explanitory
* **-a**: Analysis. Name used to determine with cc/TSelector will fun over your files. Map defined in ```./Utilities/python/SelectorTools.py``` (Name from AnalysisDatasetManager)
* **-f**: Filenames. The name of the files to be run over, in quotes seperated by commas. The filenames are those specified in the AnalysisDatasetManager, specifically with the FileInfo folder.
* **-j**: Number of Cores. Same as with Make, number of cores used to run the jobs
* **-c**: Channels: If you want to only run over a certain number of channels, put those channels in quotes seperated by commas. Default to Inclusive
* **-o**: Output File: Name of the outfile. One one is made, all the samples are put into folders after the filename
* **-b**: Hist Names: If you only want specific histograms, put names in quotes seperated by commas
* **--input_tier**: Name corresponding to the skim used on the files (ie FileInfo area)
* **--lumi**: Self Explanitory
* **--test**: Confusing naming (without context at least), but basically doesn't do Data driven backgrounds
* **--uwvv**: Legacy Ntuple format
* **--noHistConfig**: Ignore Histogram information from ADM. Generally unsafe, would ignore

So, if you want to run a basic analysis, you could write:

``` sh
./Utilities/scripts/makeHistFile.py -a Zstudy_2016 
	                                --input_tier NanoDileptonSkims 
									-s TightWithLooseVeto
									-f "dy, data_2016"
									--lumi 35.9
									--test
									-o output_zstudy.root
```
So this corresponds to running the ZSelector over the events in the dy and data_2016 files that were skimmed with NanoDileptonSkims and analyzed to TightWithLooseVeto. The Luminosity is 35.9 and we are ignoring data driven backgrounds (--test)

## Running parallel local analysis

Parallel mode currently only supports parallelizing by dataset. That is, each dataset will spawn a new process in python. This is activated by using the ```--j <num_cores>``` option to makeHistFile.py. It is not supported to parallelize locally by file, for this, using condor is suggested.

## Running anlaysis on condor

You can farm out the processing of individual files from a data set using the script [submitMakeHistFileToCondor.py](Utilities/scripts/submitMakeHistFileToCondor.py). Run ```./Utilities/scripts/submitMakeHistFileToCondor.py --help``` to see options for this script. 

Not that the condor jobs require afs in order to access some configuration files for the CMSSW/ROOT setup. Work is in progress to make this independent of afs, but currently it is expected that you use either lxplus or another machine where condor can access afs. 

NOTE: At it is not allowed to read any proxy stored below /tmp directory (which is a default directory to store the proxy). You should set this to a path on afs, for example.:

```
export X509_USER_PROXY=”/afs/cern.ch/user/${USER::1}/$USER/private/proxy/myproxy”
voms-proxy-init --voms=cms --valid=168:00
```

An example command for the ZGen analysis is

```
./Utilities/scripts/submitMakeHistFileToCondor.py -f DYm50 -a ZGen -s None --input_tier NanoAOD -d ~/work/Submit_DYtest -n 1 -q longlunch --removeUnmerged --merge ~/work/Submit_DYtest/DYtest.root 0.9 --selectorArgs theoryUnc=1 --submit
```

The name "DYm50" is specified in the AnalysisDatasetManager in FileInfo/ZGen/NanoAOD.py. The DAS path stored here is used to look up all the files for this dataset and to create a list of files passed to condor. It is expected that the files are readable with xrootd. If you create private files, they should either be published, or be in a director that is accessible with ls, e.g., eos.

### Implementing your own selector

## Running Statistical Analysis with Higgs combine

## Plotting

## Producing UWVV Ntuples

See the documentation in [UWVV](https://github.com/uwcms/UWVV]). It is possible to use other ntuples, but you will need to make extensive changes to variable names in the skimming step and branch names in the analysis code.

## Skimming UWVV Ntuples

This code is based on ["cut strings"](https://root.cern.ch/doc/v608/classTCut.html) applied to ROOT TTrees. It is driven by the script [skimNtuples.py](skimNtuples.py). This script takes as argument the analysis name and a list of selections to be applied. These selections are defined in the [Cuts folder](Cuts) of this repository. To implement a new analysis, one should create a new folder in this repository. To add a selection, a new selection.json file should be added. Follow the example of e.g. [3LooseLeptonsNoVeto.json](Cuts/WZxsec2016/3LooseLeptonsNoVeto.json). Conditions can be object specific (e.g. pt for each muon or elector) or state specific.

Run ```./skimNtuples.py --help``` for a full list of arguments to this script.

Generally you will want to farm this step out to condor, as it's real utility is to take files on hdfs, skim them, and produce output which can be copied locally. The script [farmoutNtupleSkim.py](farmoutNtupleSkim.py) is intended for this. It reads information about the input datasets from [AnalysisDatasetManager](https://github.com/kdlong/AnalysisDatasetManager) and configures jobs to be submitted to condor using ```farmoutAnalysisJobs.sh```, a script which is UW condor specific. Modifications would be required for use with other resources.

Run ```./farmoutNtupleSkim.py --help``` for a full list of arguments to this script.

An example to produce the output for the WZ inclusive analysis, with loosely IDed leptons (necessary for the fake rate step) would be

```./farmoutNtupleSkim.py -f data* -s 3MediumLeptonsNoVeto
```

This will submit jobs for each file with a name matching the pattern "data*", defined in AnalysisDatasetManager, creating a skim of events passing the [3LooseLeptonsNoVeto.json](Cuts/WZxsec2016/3LooseLeptonsNoVeto.json) selection. The script creates submit folders for each dataset, by default in the location ```/<storage>/<username>/YYYY-MM-DD_VVAnalysisJobs```, where <storage> is either /nfs_scratch or /data. It will produce output files copied to ```/store/user/<username>/VVAnalysisJobs_YYYY-MM-DD```. If you want to copy these locally, you can use the script [copyFromHdfs.py](https://github.com/kdlong/AnalysisDatasetManager/blob/master/copyFromHdfs.py).
  
``` ./copyFromHdfs.py /hdfs/store/user/<username>/...
```

## Addendum: AnalysisDatasetManager

The [AnalysisDatasetManager](https://github.com/kdlong/AnalysisDatasetManager) is a repo used for handing generic histogram and file information across softwares. Because of it's importance, it's worth explaining (may move this to the actual repo, but here for now). It is worth it to fork your own version of this repo.

This section will divided up into explaining each of the important folders in this repo and how they relate to the code

### FileInfo

The FileInfo folder contains two main parts. The first is a folder named after the Anaysis in question. As per the Z selector example above, there is a Zstudy_2016 folder for that specific study. Inside each analysis folder is a python file (all the data is stored in json-esque data format, but in python file for convience) with names corresponding to the input-tier/skim name. If our input files were skimmed by the "NanoDileptonSkim" skim for our Z study, we'd expect to see a ```./AnalysisDatasetManager/FileInfo/Zstudy_2016/NanoDileptonSkim.py ``` file that contains the information of the files. For an example of what the file should look like, [look here](https://github.com/kdlong/AnalysisDatasetManager/blob/master/FileInfo/Zstudy_2016/NanoDileptonSkim.py).

The other use is to cold monte carlo information that is used for plotting. It also allows you to put extra information about the MC's name, generator, etc so that has easy reference for later. It also allows you to put in k-factors so that doesn't have to be done by hand. All of this is locate in ```./AnalysisDatasetManager/FileInfo/montecarlo/```.

### PlotGroups

The PlotGroups gives a place to show how different MC samples will be added together in the plot. So for Zstudy_2016, there is a ```PlotGroups/Zstudy_2016.py``` that lists which files go with which overall group and assigns a name and color to it. These colors are what are specified in the ```./AnalysisDatasetManager/Styles``` folder.

### PlotObjects

Like the FileInfo directory, the PlotObjects contains folders for each analysis and a file in each analysis folder for the selection made. These files each contain information about the different plots, the initialize subcategory allowing for setting up the histogram based on the input quantites for a TH1 object and the Attributes being the ROOT functions to be run on the TH1 with it's corresponding input value. From this, the axis can be labelled and simple formatting can be done with relative ease. 

-------------------------------------------------------------------------------

With all of these understood, one can create different selection and apply different coloring convientions, have different skims for different studies, etc in one place so that one's skimming, analysis, and plotting software all have a common area to look for the configuration input.


