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
        self.correlateScaleUnc = False
        self.channels = []
        self.variations = {}
        self.normalizedVariations = []
        self.perbinVariations = {}
        self.rebin = None
        self.isMC = True
        self.isUnrolledFit = False
        self.removeZeros = True
        self.lumi = 1
        self.outputFolder = "."
        self.channelsToCombine = {}
        self.theoryVariations = {}
        self.extraCardVars = ""
        self.cardGroups = ""
        self.addOverflow = False

    def setPlotGroups(self, xsecMap):
        self.crossSectionMap = xsecMap

    def setRebin(self, rebin):
        self.rebin = rebin

    def setAddOverflow(self, overflow):
        self.addOverflow = overflow

    def setCorrelateScaleUnc(self, correlate):
        self.correlateScaleUnc = correlate

    def setRemoveZeros(self, removeZeros):
        self.removeZeros = removeZeros

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

    def setNormalizedVariations(self, variations, exclude=[]):
        self.normalizedVariations = variations
        self.setVariations(variations, exclude, True)

    def setVariations(self, variations, exclude=[], extend=False):
        if not self.processes:
            raise ValueError("No processes defined, can't set variations")
        for process in self.processes.keys():
            if process in exclude:
                self.setVariationsByProcess(process, [], extend)
                continue
            self.setVariationsByProcess(process, variations, extend)

    def getVariations(self):
        return self.variations

    def getVariationsForProcess(self, process):
        if process not in self.variations.keys():
            raise ValueError("Variations not defined for process %s" % process)
        return self.variations[process]

    def setVariationsByProcess(self, process, variations, extend=False):
        if "Up" not in variations and "Down" not in variations:
            variations = [x+y for x in variations for y in ["Up", "Down"]]
        if not extend:
            self.variations[process] = variations
        else:
            self.variations[process].extend(variations)

    def weightHistName(self, channel, process, append=""):
        fitVariable = self.getFitVariable(process)
        if append != "":
            # Removing LHE is just a hack for now!
            fitVariable = fitVariable + "_" + append
        variable = fitVariable.replace("unrolled", "2D") if self.isUnrolledFit else fitVariable
        return "_".join([variable, "lheWeights", channel])

    def setLumi(self, lumi):
        self.lumi = lumi

    def setOutputFolder(self, outputFolder):
        self.outputFolder = outputFolder
        if not os.path.isdir(outputFolder):
            os.makedirs(outputFolder)

    # Groups are muR, muF, muRmuF
    def setScaleVarGroups(self, processName, groups):
        if processName not in self.theoryVariations or 'scale' not in self.theoryVariations[processName]:
            raise ValueError("Cannot define variation groups before scale is defined for process")
        self.theoryVariations[processName]['scale']['groups'] = groups

    def addPerBinVariation(self, processName, varName, variation, correlate=True):
        if processName not in self.perbinVariations:
            self.perbinVariations[processName] = []
        self.perbinVariations[processName].append((varName, variation, correlate))

    def perBinVariationHists(self, hist, varName, var, correlation, addToCard=False):
        varUpRef = hist.Clone(hist.GetName().replace(self.fitVariable, '_'.join([self.fitVariable, varName+"Up"])))
        varDownRef = varUpRef.Clone(varUpRef.GetName().replace("Up", "Down"))

        varHists = []
        # TODO: this needs to be made more general
        if correlation and addToCard:
            self.extraCardVars += "%s    shape   1\n" % varName
        for i in range(1, hist.GetNbinsX()+1):
            if not correlation:
                varUp = varUpRef.Clone(varUpRef.GetName().replace(varName, varName+str(i)))
                varDown = varDownRef.Clone(varDownRef.GetName().replace(varName, varName+str(i)))
                self.extraCardVars += "%s%i    shape   1\n" % (varName, i)
            else:
                varUp = varUpRef
                varDown = varDownRef

            varUp.SetBinContent(i, (1.+var)*hist.GetBinContent(i))
            varUp.SetBinError(i, (1.+var)*hist.GetBinError(i))
            varDown.SetBinContent(i, 1/(1.+var)*hist.GetBinContent(i))
            varDown.SetBinError(i, 1/(1.+var)*hist.GetBinError(i))
            if i == 0 or not correlation:
                varHists.extend([varUp, varDown])
                
        return varHists

    def addTheoryVar(self, processName, varName, entries, central=0, exclude=[], specName=""):
        if varName.lower() not in ["scale", "pdf_hessian", "pdf_mc", "pdf_assymhessian", "other"]:
            raise ValueError("Invalid theory uncertainty %s. Must be type 'scale,' 'pdf,' or 'other'" % varName)
        name = "scale" if "scale" in varName.lower() else ("pdf"+specName if "pdf" in varName.lower() else specName)

        if not processName in self.theoryVariations:
            self.theoryVariations[processName] = {}

        self.theoryVariations[processName].update({ name : {
                "name" : specName,
                "entries" : entries,
                "central" : central,
                "exclude" : exclude,
                "groups" : [(3,6), (1,2), (4,8)],
                "combine" : "envelope" if varName in ["scale", "other"] else varName.replace("pdf_", ""),
            }
        })

    def addScaleBasedVar(self, processName, variation):
        if 'scale' not in self.theoryVariations[processName]:
            raise ValueError("Scale variations not defined for process %s." % processName \
                                + " Can't add additional variation." % procesName)

        scaleVars = self.theoryVariations[processName]['scale']
        if 'theoryBasedVars' not in scaleVars:
            scaleVars['theoryBasedVars'] = [variation]
        else:
            scaleVars['theoryBasedVars'].append(variation)

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
            variations.insert(0, "")
        for label, channels in self.channelsToCombine.iteritems():
            if label not in self.yields:
                self.yields[label] = {}
            for var in variations:
                name = "_".join([fitVariable, var]) if var != "" else fitVariable
                hist_name = name + "_" + channels[0]
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
        if processName in self.theoryVariations:
            plots += [self.weightHistName(c, processName) for c in self.channels]
            theoryVars = self.theoryVariations[processName]
            if 'scale' in theoryVars and 'theoryBasedVars' in theoryVars['scale']:
                thVars = theoryVars['scale']['theoryBasedVars']
                plots += [self.weightHistName(c, processName, a) for c in self.channels for a in thVars]
        return plots

    # processName needs to match a PlotGroup 
    def loadHistsForProcess(self, processName, scaleNorm=1, expandedTheory=True):
        fitVariable = self.getFitVariable(processName)
        plotsToRead = self.listOfHistsByProcess(processName, nameReplace=("unrolled", "2D"))

        group = HistTools.makeCompositeHists(self.inputFile, processName, 
                    {proc : self.crossSectionMap[proc] for proc in self.processes[processName]}, 
                    self.lumi, plotsToRead, rebin=self.rebin, overflow=False)

        if self.isUnrolledFit:
            for hist in group:
                if not "TH2" in hist.ClassName():
                    continue
                hist = HistTools.makeUnrolledHist(hist, self.unrolledBinsX, self.unrolledBinsY, overflow=self.addOverflow)
                histName = hist.GetName()
                group.Add(hist)

        #TODO:Make optional
        processedHists = []
        for chan in self.channels:
            histName = "_".join([fitVariable, chan]) if chan != "all" else fitVariable
            hist = group.FindObject(histName)
            if not hist:
                logging.warning("Failed to produce hist %s for process %s" % (histName, processName))
                continue
            if self.removeZeros and "data" not in processName.lower():
                HistTools.removeZeros(hist)
            if self.addOverflow:
                HistTools.addOverflow(hist)
            processedHists.append(histName)
            self.yields[chan].update({processName : round(hist.Integral(), 3) if hist.Integral() > 0 else 0.0001})

            if chan == self.channels[0]:
                self.yields["all"][processName] = self.yields[chan][processName]
            else:
                self.yields["all"][processName] += self.yields[chan][processName]

            if processName in self.perbinVariations:
                for varName, var, corr in self.perbinVariations[processName]:
                    map(lambda x: group.Add(x), self.perBinVariationHists(hist, varName, var, corr, chan == self.channels[0]))

            scaleHists = []
            if processName in self.theoryVariations:
                theoryVars = self.theoryVariations[processName]
                try:
                    scaleHists.extend(self.scaleHistsForProcess(group, processName, chan, expandedTheory=True))
                    if 'theoryBasedVars' in theoryVars['scale']:
                        for theoryBasedVar in theoryVars['scale']['theoryBasedVars']:
                            scaleHists.extend(self.scaleHistsForProcess(group, processName, chan, expandedTheory, theoryBasedVar))
                except ValueError as e:
                    logging.warning(e)
                    continue

                others = [k for k in theoryVars.keys() if theoryVars[k]['combine'] == 'envelope' and k != 'scale']
                for other in others:
                    group.extend(self.lheVarHistsForProcess(group, processName, chan, varName=other))

                pdfHists = []
                pdfVars = filter(lambda x: 'pdf' in x, theoryVars.keys())
                weightHist = group.FindObject(self.weightHistName(chan, processName))
                for var in pdfVars: 
                    pdfVar = theoryVars[var]
                    pdfType = "MC"
                    if "hessian" in pdfVar['combine']:
                        pdfType = "Hessian" if "assym" not in pdfVar['combine'] else "AssymHessian"

                    pdfFunction = "get%sPDFVarHists" % pdfType 
                    pdfUncScale = (1.0/1.645) if "CT18" in pdfVar['name'] else 1.0
                    args = [weightHist, pdfVar['entries'], processName, self.rebin, pdfVar['central'], pdfVar['name'], pdfUncScale]
                    if self.isUnrolledFit:
                        pdfFunction = pdfFunction.replace("get", "getTransformed3D")
                        args = args[0:1] + [HistTools.makeUnrolledHist, [self.unrolledBinsX, self.unrolledBinsY]] + args[1:]
                    updatePdfs = getattr(HistTools, pdfFunction)(*args)
                    pdfHists += updatePdfs

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
                        theoryVarLabels.extend(re.findall("_".join([fitvar, "(.*pdf.*)", chan]), h.GetName()))
                        theoryVarLabels.extend(re.findall("_".join([fitvar, "(.*QCDscale.*)", chan]), h.GetName()))
                    self.variations[processName].extend(set(theoryVarLabels))

                normVarNames = []
                map(lambda x: normVarNames.extend([x+"Up", x+"Down"]), self.normalizedVariations)
                for normvar in normVarNames:
                    varHistName = histName.replace(chan, "_".join([normvar, chan]))
                    hist = group.FindObject(varHistName)
                    if not hist:
                        logging.warning("Failed to find hist %s in group %s. Skipping" % (varHistName, group.GetName()))
                        continue
                    hist.Scale(self.yields[chan][processName]/hist.Integral())

        if self.addOverflow:
            map(HistTools.addOverflow, filter(lambda x: (x.GetName() not in processedHists), group))
        if self.removeZeros and "data" not in group.GetName().lower():
            map(HistTools.removeZeros, filter(lambda x: (x.GetName() not in processedHists), group))
        #TODO: You may want to combine channels before removing zeros
        if self.channelsToCombine.keys():
            self.combineChannels(group, processName)

        self.histData[processName] = group
        
    def lheVarHistsForProcess(self, group, processName, chan, varName):
        weighthist_name = self.weightHistName(chan, processName)
        weightHist = group.FindObject(weighthist_name)
        if not weightHist:
            raise ValueError("Failed to find %s. Skipping" % weighthist_name)
        var = self.theoryVariations[processName][varName]
        if not self.isUnrolledFit:
            hists,name = HistTools.getLHEWeightHists(weightHist, var['entries'], "", varName)
        else:
            hists,name = HistTools.getTransformed3DLHEHists(weightHist, HistTools.makeUnrolledHist,
                [self.unrolledBinsX, self.unrolledBinsY], var['entries'], "", varName)
        
        # In this case, should be central, varCentral, varUp, varDown
        if len(var['entries']) == 4:
            for i in range(hists[0].GetNbinsX()+1):
                central_ratio = hists[0].GetBinContent(i)/hists[1].GetBinContent(i) if hists[1].GetBinContent(i) > 0 else 1
                hists[2].SetBinContent(i, hists[2].GetBinContent(i)*central_ratio)
                hists[3].SetBinContent(i, hists[3].GetBinContent(i)*central_ratio)
            hists = [hists[2], hists[3]]

        hists[0] = HistTools.rebinHist(hists[0], name, self.rebin)
        hists[1] = HistTools.rebinHist(hists[1], name.replace("Up", "Down"), self.rebin)
        return hists

    def scaleHistsForProcess(self, group, processName, chan, expandedTheory, append=""):
        weighthist_name = self.weightHistName(chan, processName)
        weightHist = group.FindObject(weighthist_name)
        if not weightHist:
            raise ValueError("Failed to find %s. Skipping" % weighthist_name)
        if 'scale' not in self.theoryVariations[processName]:
            return []
        scaleVars = self.theoryVariations[processName]['scale']
        
        scaleHists = HistTools.getScaleHists(weightHist, processName, self.rebin, 
                entries=scaleVars['entries'], 
                exclude=scaleVars['exclude'], 
                central=scaleVars['central']) if not self.isUnrolledFit else \
            HistTools.getTransformed3DScaleHists(weightHist, HistTools.makeUnrolledHist,
                    [self.unrolledBinsX, self.unrolledBinsY], processName,
                entries=scaleVars['entries'], 
                exclude=scaleVars['exclude'])

        if expandedTheory:
            name = processName if not self.correlateScaleUnc else ""
            expandedScaleHists = HistTools.getExpandedScaleHists(weightHist, name, self.rebin, 
                    entries=scaleVars['entries'], 
                    pairs=scaleVars['groups'], 
                ) if not self.isUnrolledFit else \
                HistTools.getTransformed3DExpandedScaleHists(weightHist, 
                        HistTools.makeUnrolledHist,
                    [self.unrolledBinsX, self.unrolledBinsY], name,
                    entries=scaleVars['entries'], 
                    pairs=scaleVars['groups'], 
                )
                
            scaleHists.extend(expandedScaleHists)
        return scaleHists

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
        chan_dict = self.yields[chan].copy()
        chan_dict.update(extraArgs)
        for key, value in extraArgs.iteritems():
            if "yield:" in value:
                chan_dict[key] = chan_dict[value.replace("yield:", "")]
        chan_dict["nuisances"] = nuisances
        chan_dict["fit_variable"] = self.fitVariable
        chan_dict["output_file"] = self.outputFile.GetName()
        chan_dict["card_append"] = self.extraCardVars + "\n\n" + self.cardGroups 
        outputCard = self.templateName.split("/")[-1].format(channel=chan, label=label) 
        outputCard = outputCard.replace("template", outlabel)
        outputCard = outputCard.replace("__", "_")
        ConfigureJobs.fillTemplatedFile(self.templateName.format(channel=chan, label=label),
            "/".join([self.outputFolder, outputCard]),
            chan_dict
        )
        chan_dict.pop("card_append")
        print chan_dict

