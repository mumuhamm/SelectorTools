import logging
import re
import ConfigureJobs
import HistTools
import OutputTools
from prettytable import PrettyTable
import os
import ROOT

class CombineCardTools(object):
    def __init__(self):
        self.fitVariable = ""
        self.fitVariableAppend = {}
        self.processes = []
        self.yields = {}
        self.histData = {}
        self.crossSectionMap = {}
        self.outputFile = ""
        self.templateName = ""
        self.channels = []
        self.variations = {}
        self.rebin = None
        self.isMC = True
        self.isUnrolledFit = False
        self.lumi = 1
        self.outputFolder = "."
        self.channelsToCombine = {}
        self.theoryVariations = {}

    def setPlotGroups(self, xsecMap):
        self.crossSectionMap = xsecMap

    def setRebin(self, rebin):
        self.rebin = rebin

    def setUnrolled(self, binsx, binsy):
        self.isUnrolledFit = True
        self.unrolledBinsX = binsx
        self.unrolledBinsY = binsy

    def setCrosSectionMap(self, xsecMap):
        self.crossSectionMap = xsecMap

    # Map of plot groups and members (individual processes)
    def setProcesses(self, processes):
        self.processes = processes

    def setFitVariableAppend(self, process, append):
        self.fitVariableAppend[process] = append

    def setFitVariable(self, variable):
        self.fitVariable = variable

    def setVariations(self, variations, exclude=[]):
        if not self.processes:
            raise ValueError("No processes defined, can't set variations")
        for process in self.processes.keys():
            if process in exclude:
                self.setVariationsByProcess(process, [])
                continue
            self.setVariationsByProcess(process, variations)

    def getVariations(self):
        return self.variations

    def getVariationsForProcess(self, process):
        if process not in self.variations.keys():
            raise ValueError("Variations not defined for process %s" % process)
        return self.variations[process]

    def setVariationsByProcess(self, process, variations):
        if "Up" not in variations and "Down" not in variations:
            variations = [x+y for x in variations for y in ["Up", "Down"]]
        self.variations[process] = variations

    def weightHistName(self, channel, process):
        fitVariable = self.getFitVariable(process)
        variable = fitVariable.replace("unrolled", "2D") if self.isUnrolledFit else fitVariable
        return "_".join([variable, "lheWeights", channel])

    def setLumi(self, lumi):
        self.lumi = lumi

    def setOutputFolder(self, outputFolder):
        self.outputFolder = outputFolder
        if not os.path.isdir(outputFolder):
            os.makedirs(outputFolder)

    def setScaleVarGroups(self, processName, groups):
        if processName not in self.theoryVariations or 'scale' not in self.theoryVariations[processName]:
            raise ValueError("Cannot define variation groups before scale is defined for process")
        self.theoryVariations[processName]['scale']['groups'] = groups


    def addTheoryVar(self, processName, varName, entries, central=0, exclude=[], specName=""):
        if "scale" not in varName.lower() and "pdf" not in varName.lower():
            raise ValueError("Invalid theory uncertainty %s. Must be type 'scale' or 'pdf'" % varName)
        name = "scale" if "scale" in varName.lower() else ("pdf"+specName)

        if not processName in self.theoryVariations:
            self.theoryVariations[processName] = {}

        self.theoryVariations[processName].update({ name : {
                "name" : specName,
                "entries" : entries,
                "central" : central,
                "exclude" : exclude,
                "groups" : [(3,6), (1,2), (4,8)],
                "combine" : "envelope" if name == "scale" else varName.replace("pdf_", ""),
            }
        })

    def getRootFile(self, rtfile, mode=None):
        if type(rtfile) == str:
            if mode:
                return ROOT.TFile.Open(rtfile, mode)
            else:
                return ROOT.TFile.Open(rtfile)
        return rtfile

    def setTemplateFileName(self, templateName):
        self.templateName = templateName

    def setOutputFile(self, outputFile):
        self.outputFile = self.getRootFile("/".join([self.outputFolder, outputFile]), "RECREATE")

    def setInputFile(self, inputFile):
        self.inputFile = self.getRootFile(inputFile)

    def setChannels(self, channels):
        self.channels = channels
        for chan in self.channels + ["all"]:
            self.yields[chan] = {}

    def processHists(self, processName):
        return self.histData[processName] 

    def getFitVariable(self, process, replace=None):
        fitvar = self.fitVariable.replace(*replace) if replace else self.fitVariable
        if process not in self.fitVariableAppend:
            return fitvar
        return "_".join([fitvar, self.fitVariableAppend[process]])

    def setCombineChannels(self, groups):
        self.channelsToCombine = groups

    def combineChannels(self, group, processName, central=True):
        variations = self.variations[group.GetName()][:]
        fitVariable = self.getFitVariable(group.GetName())
        if central:
            variations.append("")
        for label, channels in self.channelsToCombine.iteritems():
            if label not in self.yields:
                self.yields[label] = {}
            for var in variations:
                name = "_".join([fitVariable, var]) if var != "" else fitVariable
                hist_name = name + "_" + channels[0]
                tmphist = group.FindObject(hist_name)
                if not tmphist:
                    logging.warning("Failed to find hist %s in group %s. Skipping" % (hist_name, group.GetName()))
                    continue
                hist = tmphist.Clone("_".join([name, label]))
                del tmphist
                ROOT.SetOwnership(hist, False)
                group.Add(hist) 
                for chan in channels[1:]:
                    chan_hist = group.FindObject(name + "_" + chan)
                    hist.Add(chan_hist)
                if var == "":
                    self.yields[label][processName] = round(hist.Integral(), 3)

    def listOfHistsByProcess(self, processName, nameReplace=None):
        if self.fitVariable == "":
            raise ValueError("Must declare fit variable before defining plots")
        fitVariable = self.getFitVariable(processName, replace=nameReplace)
        plots = ["_".join([fitVariable, chan]) for chan in self.channels]
        variations = self.getVariationsForProcess(processName)
        plots += ["_".join([fitVariable, var, c]) for var in variations for c in self.channels]
        if processName in self.theoryVariations.keys():
            plots += [self.weightHistName(c, processName) for c in self.channels]
        return plots

    # processName needs to match a PlotGroup 
    def loadHistsForProcess(self, processName, scaleNorm=1, expandedTheory=True):
        plotsToRead = self.listOfHistsByProcess(processName, nameReplace=("unrolled", "2D"))

        group = HistTools.makeCompositeHists(self.inputFile, processName, 
                    {proc : self.crossSectionMap[proc] for proc in self.processes[processName]}, 
                    self.lumi, plotsToRead, rebin=self.rebin, overflow=False)

        if self.isUnrolledFit:
            for hist in group:
                if not plot.InheritsFrom("TH2"):
                    continue
                hist = HistTools.makeUnrolledHist(plot, self.unrolledBinsX, self.unrolledBinsY)
                histName = hist.GetName()
                group.Add(hist)

        fitVariable = self.getFitVariable(processName)
        #TODO:Make optional
        processedHists = []
        for chan in self.channels:
            histName = "_".join([fitVariable, chan]) if chan != "all" else fitVariable
            hist = group.FindObject(histName)
            if not hist:
                logging.warning("Failed to produce hist %s for process %s" % (histName, processName))
                continue
            #TODO: Make optional
            if "data" not in processName.lower():
                HistTools.removeZeros(hist)
            HistTools.addOverflow(hist)
            processedHists.append(histName)
            self.yields[chan].update({processName : round(hist.Integral(), 3) if hist.Integral() > 0 else 0.0001})
            print self.yields

            if chan == self.channels[0]:
                self.yields["all"][processName] = self.yields[chan][processName]
            else:
                self.yields["all"][processName] += self.yields[chan][processName]

            if processName in self.theoryVariations:
                weightHist = group.FindObject(self.weightHistName(chan, processName))
                if not weightHist:
                    logging.warning("Failed to find %s. Skipping" % self.weightHistName(chan, processName))
                    continue
                theoryVars = self.theoryVariations[processName]
                scaleHists = HistTools.getScaleHists(weightHist, processName, self.rebin, 
                        entries=theoryVars['scale']['entries'], 
                        exclude=theoryVars['scale']['exclude'], 
                        central=theoryVars['scale']['central']) if not self.isUnrolledFit else \
                    HistTools.getTransformed3DScaleHists(weightHist, HistTools.makeUnrolledHist,
                            [self.unrolledBinsX, self.unrolledBinsY], processName,
                        entries=theoryVars['scale']['entries'], 
                        exclude=theoryVars['scale']['exclude'])
                if expandedTheory:
                    expandedScaleHists = HistTools.getExpandedScaleHists(weightHist, processName, self.rebin, 
                            entries=theoryVars['scale']['entries'], 
                            pairs=theoryVars['scale']['groups'], 
                        ) if not self.isUnrolledFit else \
                        HistTools.getTransformed3DExpandedScaleHists(weightHist, 
                                HistTools.makeUnrolledHist,
                            [self.unrolledBinsX, self.unrolledBinsY], processName,
                            entries=theoryVars['scale']['entries'], 
                            pairs=theoryVars['scale']['groups'], 
                        )
                        
                    scaleHists.extend(expandedScaleHists)

                pdfHists = []
                pdfVars = filter(lambda x: 'pdf' in x, theoryVars.keys())
                for var in pdfVars: 
                    pdfVar = theoryVars[var]
                    pdfType = "MC"
                    if "hessian" in pdfVar['combine']:
                        pdfType = "Hessian" if "assym" not in pdfVar['combine'] else "AssymHessian"

                    pdfFunction = "get%sPDFVarHists" % pdfType 
                    args = [weightHist, pdfVar['entries'], processName, self.rebin, pdfVar['central'], pdfVar['name']]
                    if self.isUnrolledFit:
                        pdfFunction = pdfFunction.replace("get", "getTransformed3D")
                        args = args[0:1] + [HistTools.makeUnrolledHist, [self.unrolledBinsX, self.unrolledBinsY]] + args[1:]

                    pdfHists += getattr(HistTools, pdfFunction)(*args)

                    if expandedTheory and pdfVar['name'] == 'NNPDF31':
                        args.pop(len(args)-1)
                        pdfFunction = HistTools.getAllSymHessianHists if not self.isUnrolledFit else HistTools.getTransformed3DAllSymHessianHists
                        allPdfHists = pdfFunction(*args)
                        pdfHists.extend(allPdfHists)

                        if not self.isUnrolledFit:
                            cenHist, _ = HistTools.getLHEWeightHists(weightHist, pdfVar['entries'][:1], processName, "pdf", self.rebin)
                        else:
                            cenHist = HistTools.getAllTransformed3DHists(weightHist, HistTools.makeUnrolledHist, 
                                    [self.unrolledBinsX, self.unrolledBinsY], processName, pdfVar['entries'][:1])
                        if len(cenHist) and cenHist[0]:
                            pdfHists.append(cenHist[0].Clone())
                group.extend(scaleHists+pdfHists)

                if chan == self.channels[0]:
                    theoryVarLabels = []
                    for h in group:
                        fitvar = self.getFitVariable(processName)
                        theoryVarLabels.extend(re.findall("_".join([fitvar, "(.*)", chan]), h.GetName()))
                    self.variations[processName].extend(set(theoryVarLabels))

        #TODO: Make optional
        map(HistTools.addOverflow, filter(lambda x: (x.GetName() not in processedHists), group))
        if "data" not in group.GetName().lower():
            map(HistTools.removeZeros, filter(lambda x: (x.GetName() not in processedHists), group))
        #TODO: You may want to combine channels before removing zeros
        self.combineChannels(group, processName)

        self.histData[processName] = group

    #def makeNuisanceShapeOnly(process, central_hist, uncertainty, channels):
    #    for chan in channels:
    #        central_hist = self.histData[process].FindObject("_".join([central_hist, chan]))
    #        var_hist = self.histData[process].FindObject("_".join([central_hist, uncertainty, chan]))

    # It's best to call this function for process, otherwise you can end up
    # storing many large histograms in memory
    def writeProcessHistsToOutput(self, processName):
        if processName not in self.histData.keys() or not self.histData[processName]:
            raise ValueError("Hists for process %s not found" % processName)
        processHists = self.histData[processName]
        OutputTools.writeOutputListItem(processHists, self.outputFile)
        processHists.Delete()
        
    def writeMetaInfo(self):
        OutputTools.addMetaInfo(self.outputFile)

    def writeCards(self, chan, nuisances, label="", outlabel="", extraArgs={}):
        print self.yields
        chan_dict = self.yields[chan].copy()
        print chan_dict
        chan_dict.update(extraArgs)
        for key, value in extraArgs.iteritems():
            if "yield:" in value:
                chan_dict[key] = chan_dict[value.replace("yield:", "")]
        chan_dict["nuisances"] = nuisances
        chan_dict["fit_variable"] = self.fitVariable
        chan_dict["output_file"] = self.outputFile.GetName()
        outputCard = self.templateName.split("/")[-1].format(channel=chan, label=label) 
        outputCard = outputCard.replace("template", outlabel)
        outputCard = outputCard.replace("__", "_")
        ConfigureJobs.fillTemplatedFile(self.templateName.format(channel=chan, label=label),
            "/".join([self.outputFolder, outputCard]),
            chan_dict
        )

