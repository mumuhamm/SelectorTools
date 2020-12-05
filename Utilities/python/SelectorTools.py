#!/usr/bin/env python
import ROOT
import glob
import datetime
import ConfigureJobs, OutputTools
import sys
import os
import multiprocessing
import subprocess
import logging
import numpy

class SelectorDriver(object):
    def __init__(self, analysis, selection, input_tier, year):
        selector_map = {
            "WZxsec2016" : "WZSelector",
            "Zstudy" : "ZSelector",
            "LowPileupZ" : "LowPileupZSelector",
            "LowPileupW" : "LowPileupWSelector",
            "Zstudy_2016" : "ZSelector",
            "Zstudy_2017" : "ZSelector",
            "ZZGen" : "ZZGenSelector",
            "WGen" : "WGenSelector",
            "ZGen" : "ZGenSelector",
            "TTT" : "TTTSelector",
            "ThreeLep" : "ThreeLepSelector",
            "Eff" : "Efficiency",
            "Efficiency" : "Efficiency",
        }

        self.select = None
        self.subanalysis = None
        if analysis.find(":") != -1:
            analysis, self.subanalysis = analysis.split(':')
        self.analysis = analysis
        self.selection = selection
        self.input_tier = input_tier
        if analysis not in selector_map.keys():
            raise ValueError("Analysis does not point to " \
                "a defined selector. Please edit "
                "Utilities/python/SelectorTools.py to add it.")
        self.selector_name = selector_map[analysis]
        self.addSumweights = True
        self.ntupleType = "NanoAOD"
        self.year = year
        self.numCores = 1
        self.channels = ["Inclusive"]
        self.outfile_name = "temp.root"
        self.outfile = None
        self.datasets = {}
        self.regions = {}
        self.maxEntries = -1
        self.maxFiles = -1
        self.nProcessed = 0
        self.compress = False

    # Needed to parallelize class member function, see
    # https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map
    def __call__(self, args):
        self.select = getattr(ROOT, self.selector_name)()
        self.processDatasetHelper(args)

    def tempfileName(self):
        return "temp_%s_%s" % (multiprocessing.current_process().name, self.outfile_name.split("/")[-1])

    def setChannels(self, channels):
        self.channels = channels

    def setMaxEntries(self, maxEntries):
        self.maxEntries = maxEntries

    def setMaxFiles(self, maxFiles):
        self.maxFiles = maxFiles

    def setAddSumWeights(self, addSumWeights):
        self.addSumweights = addSumWeights

    def outputFile(self):
        return self.outfile
    
    def setCompress(self, compress):
        self.compress = compress

    def updateCurrentFile(self, currfile_name):
        self.current_file = ROOT.TFile.Open(currfile_name, "update")
        if self.compress:
            self.current_file.SetCompressionAlgorithm(2)
            self.current_file.SetCompressionLevel(5)
    
    def setOutputfile(self, outfile_name):
        if self.outfile:
            self.outfile.Close()
        self.outfile_name = outfile_name
        self.outfile = ROOT.gROOT.FindObject(outfile_name)
        if not self.outfile:
            self.outfile = ROOT.TFile.Open(outfile_name, "recreate")
        if self.compress:
            self.outfile.SetCompressionAlgorithm(2)
            self.outfile.SetCompressionLevel(5)

        self.current_file = self.outfile

    def addTNamed(self, name, title):
        if not self.inputs:
            raise ValueError("Input list is not defined!")
        named = self.inputs.FindObject(name)
        if named:
            named.SetTitle(title)
        else:
            named = ROOT.TNamed(name, title)
            self.inputs.Add(named)

    def setInputs(self, inputs):
        self.inputs = ROOT.TList()
        for inp in inputs:
            self.inputs.Add(inp)
        self.addTNamed("ntupleType", self.ntupleType)
        self.addTNamed("selection", self.selection)
        self.addTNamed("year", self.year)

    def setSelection(self, selection):
        self.selection = selection
        
    def setInputTier(self, input_tier):
        self.input_tier = input_tier
        
    def setNtupeType(self, ntupleType):
        self.ntupleType = ntupleType
        self.addTNamed("ntupleType", self.ntupleType)

    def setNumCores(self, numCores):
        self.numCores = numCores

    def setFileList(self, list_of_files, nPerJob, jobNum):
        if not os.path.isfile(list_of_files):
            raise ValueError("%s is not a valid file." % list_of_files)
        filelist = [f.split("#")[0].strip() for f in open(list_of_files).readlines()]
        # Remove empty/commented lines
        filelist = filter(lambda  x: len(x) > 2, filelist)
        nPerJob = int(nPerJob)
        if nPerJob < 1:
            raise ValueError("Number of files per job must be >= 1.")
        jobNum = int(jobNum)
        maxNum = len(filelist)
        firstEntry = nPerJob*jobNum
        if firstEntry > maxNum:
            raise ValueError("The first file to process (nPerJob*jobNum) = (%i*%i)" % (nPerJob, jobNum) \
                    + " is greater than the number of entries in file %s (%s)." % (list_of_files, maxNum))
        lastEntry = min(nPerJob*(jobNum+1), maxNum)
        
        for line in filelist[firstEntry:lastEntry]:
            if "@" not in line:
                dataset = "Unknown"
                file_path = line
            else:
                # Intended for running specified files, use the format name:file
                dataset, file_path = line.split("@")
            if dataset not in self.datasets.keys():
                self.datasets[dataset] = [file_path]
            else:
                self.datasets[dataset].append(file_path)

    def setDatasetRegions(self, regions):
        regionSets = [i.strip() for i in regions.split(";")]
        for region in regionSets:
            process, regions = [i.strip() for i in region.split("=")]
            label = process
            if "__" in label:
                label, tag = process.split("__")
                tag = "__" + tag
            self.regions[process] = ["_".join([label, i.strip()+tag]) for i in regions.split(",")]
    
    def unsetDatasetRegions(self):
        self.regions = {}

    # If you want a completely new set of data, call this. All duplicates will be overwritten,
    # but if you don't define a new dataset properly, it would be read from the old if not cleared
    def clearDatasets(self):
        self.datasets = {}

    def setDatasets(self, datalist):
        analysis = self.subanalysis if self.subanalysis else self.analysis
        datasets = ConfigureJobs.getListOfFiles(datalist, self.input_tier, analysis=analysis)
        
        for dataset in datasets:
            if "@" in dataset:
                dataset, file_path = [f.strip() for f in dataset.split("@")]
            elif ".root" in dataset[-5:]:
                file_path = dataset
                dataset = "Unknown"
            else:
                try:
                    file_path = ConfigureJobs.getInputFilesPath(dataset, self.input_tier, analysis)
                except ValueError as e:
                    logging.warning(e)
                    continue

            self.datasets[dataset] = [file_path]

    def expandDatasetFilePaths(self, nsplits):
        nFiles = 0 
        for dataset, file_path in self.datasets.iteritems():
            maxPerSet = self.maxFiles/len(self.datasets) 
            if dataset == self.datasets.keys()[-1]:
                maxPerSet = self.maxFiles-nFiles
            files = []
            for f in file_path:
                files.extend([f] if os.path.isfile(f) else glob.glob(f))
            self.datasets[dataset] = numpy.array_split(files[:self.maxFiles] if self.maxFiles > 0 else files, nsplits)
            nFiles += len(self.datasets[dataset])

    def applySelector(self):
        for chan in self.channels:
            self.addTNamed("channel", chan)
            logging.info("Processing channel %s" % chan)
            if self.numCores > 1:
                self.processParallelByDataset(self.datasets, chan)
            else: 
                for dataset, file_path in self.datasets.iteritems():
                    self.processDataset(dataset, file_path, chan)
        if len(self.channels) > 1 and self.numCores > 1:
            tempfiles = [self.outfile_name.replace(".root", "_%s.root" % c) for c in self.channels]
            self.combineParallelFiles(tempfiles, "Inclusive")

    def isBackground(self):
        self.selector_name = self.selector_name.replace("Selector", "BackgroundSelector")

    def processDataset(self, dataset, file_path, chan):
        self.select = getattr(ROOT, self.selector_name)()
        logging.info("Processing dataset %s" % dataset)
        self.select.SetInputList(self.inputs)
        self.addTNamed("name", dataset)
        if dataset in self.regions:
            vec = ROOT.std.vector("string")()
            for i in self.regions[dataset]:
                vec.push_back(i)
            self.select.addSubprocesses(vec)
        # Only add for one channel
        addSumweights = self.addSumweights and self.channels.index(chan) == 0 and "data" not in dataset
        if addSumweights:
            # Avoid accidentally combining sumweights across datasets
            currfile_name = self.current_file.GetName()
            self.current_file.Close()
            sumweights_hist = ROOT.gROOT.FindObject("sumweights")
            if sumweights_hist:
                sumweights_hist.Delete()
            if not self.outfile:
                self.outfile = ROOT.TFile.Open(self.outfile_name)
            sumweights_hist = self.outfile.Get("%s/sumweights" % dataset)
            
            if not sumweights_hist:
                sumweights_hist = ROOT.TH1D("sumweights", "sumweights", 1000, 0, 1000)
            sumweights_hist.SetDirectory(ROOT.gROOT)
            if currfile_name:
                self.updateCurrentFile(currfile_name)
        self.processLocalFiles(file_path, addSumweights, chan)
        self.collectOutput(dataset, chan)

    def collectOutput(self, dataset, chan):
        addSumweights = self.addSumweights and self.channels.index(chan) == 0 and "data" not in dataset
        output_list = self.select.GetOutputList()
        processes = [dataset] + (self.regions[dataset] if dataset in self.regions else [])
        self.writeOutput(output_list, chan, processes, dataset, addSumweights)

        if self.current_file != self.outfile:
            self.current_file.Close()
        return True

    def writeOutput(self, output_list, chan, processes, dataset, addSumweights):
        sumweights_hist = ROOT.gROOT.FindObject("sumweights")

        # The file closing messes up the sumweights when its taken directly from the file
        if self.numCores > 1:
            self.outfile.Close()
            chanNum = self.channels.index(chan)
            self.updateCurrentFile(self.tempfileName())
        if not self.current_file:
            self.current_file = ROOT.TFile.Open(self.outfile_name)

        for process in processes:
            dataset_list = output_list.FindObject(process)
            if not dataset_list or dataset_list.ClassName() != "TList":
                logging.warning("No output found for process %s of dataset %s" % (process, dataset))
                dataset_list = output_list.FindObject("Unknown") if process == dataset else None
                if dataset_list and dataset_list.ClassName() == "TList":
                    logging.warning('Falling back to dataset "Unknown"')
                else:
                    logging.warning('Skipping process %s for dataset %s' % (process, dataset))
                    return False
            if addSumweights:
                if not sumweights_hist:
                    logging.warning("Failed to find sumweights for process %s" % process)
                else:
                    dataset_list.Add(sumweights_hist.Clone())
            OutputTools.writeOutputListItem(dataset_list, self.current_file)
            map(lambda x: x.Delete(), dataset_list)
            del dataset_list
        del output_list

    def getFileNames(self, file_path):
        xrootd = "store" in file_path.split("/")[0:4]
        xrootd_user = "/store/user" in file_path.replace("/eos/cms", "")[:7]
        if not (xrootd or os.path.isfile(file_path) or len(glob.glob(file_path.rsplit("/", 1)[0]))):
            raise ValueError("Invalid path! Skipping dataset. Path was %s" 
                % file_path)

        if (xrootd and not xrootd_user):
            xrd = 'root://%s/' % ConfigureJobs.getXrdRedirector(file_path)
            allfiles = glob.glob(file_path) 
            if not allfiles:
                return [xrd + file_path]
            else:
                xrd_path = lambda x: "".join([xrd, "/store", x.rsplit("/store")[1]])
                return [xrd_path(f) for f in allfiles]

        filenames =  glob.glob(file_path) if not xrootd_user else \
                ConfigureJobs.getListOfHDFSFiles(file_path)
        filenames = ['root://cmsxrootd.hep.wisc.edu/' + f if "/store/user" in f[0:12] else f for f in filenames]
        return filenames

    def getTreeName(self, chan):
        # TODO: Fix this! This is an extremely ineffient way to separate the eemm and mmee
        # since it involves reading the file an extra time
        channel = chan if "mmee" not in chan else chan.replace("mmee", "eemm")
        return ("%s/ntuple" % channel) if self.ntupleType == "UWVV" else "Events"

    def combineParallelFiles(self, tempfiles, chan):
        tempfiles = filter(os.path.isfile, tempfiles)
        outfile = self.outfile_name
        if chan != "Inclusive":
            outfile = self.outfile_name.replace(".root", "_%s.root" % chan)
        rval = subprocess.call(["hadd", "-k", "-f", "-j", str(self.numCores), outfile] + tempfiles)
        if rval == 0:
            map(os.remove, tempfiles)
        else:
            raise RuntimeError("Failed to collect data from parallel run")

    def processParallelByDataset(self, datasets, chan):
        self.expandDatasetFilePaths(self.numCores)
        expanded_datasets = [[d, f, chan] for d, files in datasets.iteritems() for f in files]
        logging.debug(expanded_datasets)
        p = multiprocessing.Pool(processes=self.numCores)
        tempfiles = glob.glob(self.tempfileName().replace("MainProcess", "PoolWorker*"))
        map(os.remove, tempfiles)
        p.map(self, expanded_datasets)
        # Store arrays in temp files, since it can get way too big to keep around in memory
        tempfiles = glob.glob(self.tempfileName().replace("MainProcess", "PoolWorker*"))
        p.close()
        self.combineParallelFiles(tempfiles, chan)

    # Pool.map can only take in one argument, so expand the array
    def processDatasetHelper(self, args):
        self.processDataset(*args)

    def processLocalFiles(self, file_path, addSumweights, chan,):
        filenames = []
        for entry in file_path:
            filenames.extend(self.getFileNames(entry))

        for i, filename in enumerate(filenames):
            self.nProcessed += self.processFile(filename, addSumweights, chan, i+1)
            logging.info("Processed %i events after %i files" % (self.nProcessed, i+1))
            if (self.maxEntries > 0 and self.nProcessed >= self.maxEntries) or (self.maxFiles > 0 and i > self.maxFiles):
                break
                
    def processFile(self, filename, addSumweights, chan, filenum=1):
        rtfile = ROOT.TFile.Open(filename)
        if not rtfile or not rtfile.IsOpen() or rtfile.IsZombie():
            raise IOError("Failed to open file %s!" % filename)
        tree_name = self.getTreeName(chan)
        tree = rtfile.Get(tree_name)
        if not tree:
            raise ValueError(("tree %s not found for file %s. " \
                    "Either the file is corrupted or the ntupleType (%s) is wrong.") 
                % (tree_name, filename, self.ntupleType)
            )
        logging.debug("Processing tree %s for file %s." % (tree.GetName(), rtfile.GetName()))
        toprocess = self.maxEntries-self.nProcessed
        if toprocess > 0:
            tree.Process(self.select, "", toprocess)
        else:
            tree.Process(self.select, "")
        logging.debug("Processed file %s with selector %s." % (filename, self.select.GetName()))
        if addSumweights:
            self.fillSumweightsHist(rtfile, filenum)
            logging.debug("Added sumweights hist.")
        entries = min(tree.GetEntries(), toprocess) if toprocess > 0 else tree.GetEntries()
        rtfile.Close()
        return entries

    # You can use filenum to index the files and sum separately, but it's not necessary
    def fillSumweightsHist(self, rtfile, filenum=1):
        sumWeightsType = "fromTree"
        weightSignOnly = filter(lambda x: "wSignOnly" in x.GetName(), self.inputs)
        wSuppress = filter(lambda x: "wSuppress" in x.GetName(), self.inputs)
        weightSignOnly = weightSignOnly[0].GetVal() if weightSignOnly else False
        wSuppress = wSuppress[0].GetVal() if wSuppress else 0

        if self.ntupleType == "NanoAOD":
            nevents_branch = "genEventCount"
            meta_tree_name = "Runs"
            sumweights_branch = "genEventSumw" 
            if not rtfile.Get(meta_tree_name).GetBranch(sumweights_branch):
                sumweights_branch += "_"
                nevents_branch += "_"
            if weightSignOnly or wSuppress:
                meta_tree_name = "Events"
                sumweights_branch = "genWeight" 
                nevents_branch = "Entries$"
        elif self.ntupleType == "Bacon":
            sumWeightsType = "fromHist"
            weightshist_name = "hGenWeights"
        else:
            nevents_branch = ""
            sumweights_branch = "summedWeights"
            meta_tree_name = "metaInfo/metaInfo"

        ROOT.gROOT.cd()
        sumweights_hist = ROOT.gROOT.FindObject("sumweights")
        if sumWeightsType == "fromTree":
            meta_tree = rtfile.Get(meta_tree_name)
            tmplabel = sumweights_hist.GetName()+"_i"
            tmpweights_hist = sumweights_hist.Clone(tmplabel)

            draw_weight = sumweights_branch 
            if weightSignOnly:
                draw_weight = "(%s > 0 ? 1 : -1)" % draw_weight
            if wSuppress:
                draw_weight = "%s*(%s < %s)" % (draw_weight, sumweights_branch, wSuppress)
            if self.maxEntries and self.maxEntries > 0:
                logging.warning("maxEntries is a subset of all the events in the file." \
                    " The sumweights hist will be scaled to reflect this, but this is NOT exact!")
                draw_weight += "*%i/%s" % (self.maxEntries, nevents_branch)

            meta_tree.Draw("%i>>%s" % (filenum, tmplabel), draw_weight)
            sumweights_hist.Add(tmpweights_hist)
        elif sumWeightsType == "fromHist":
            new_sumweights_hist = rtfile.Get(weightshist_name)
            if sumweights_hist and new_sumweights_hist:
                sumweights_hist.Add(new_sumweights_hist)
                sumweights_hist.SetDirectory(ROOT.gROOT)
                ROOT.SetOwnership(sumweights_hist, False)
        logging.debug("Sumweights is %0.2f" % sumweights_hist.Integral())
