import numpy
import array
import ROOT
import logging
import math
import numpy

def getDifference(fOut, name, dir1, dir2, ratioFunc=None):
    differences = ROOT.TList()
    differences.SetName(name)
    for histname in [i.GetName() for i in fOut.Get(dir1).GetListOfKeys()]:
        if histname == "sumweights": continue
        hist1 = fOut.Get("/".join([dir1, histname]))
        hist2 = fOut.Get("/".join([dir2, histname]))
        if hist1 and hist2:
            diff = hist1.Clone()
            diff.Add(hist2, -1)
        elif not hist1:
            logging.warning("Hist %s was not produced for " \
                "dataset(s) %s" % (histname, dir1))
        elif not hist2:
            logging.warning("Hist %s was not produced for " \
                "dataset(s) %s" % (histname, dir2))
        differences.Add(diff)
    if ratioFunc:
        ratios = ratioFunc(differences)
        for ratio in ratios:
            differences.Add(ratio) 
    return differences

def makeUnrolledHist(init_2D_hist, xbins, ybins, name=""):
    nbins = (len(xbins)-1)*(len(ybins)-1)
    hists_half_rolled = []
    for i in range(len(ybins)-1):
        ylower = ybins[i]
        yupper = ybins[i+1]
        hist_name = ("%s_%0.1fTo%0.1f" % (init_2D_hist.GetName(), ylower, yupper)).replace(".","p")
        # Global bin number = nbinx + nbinsx*nbiny
        lower_bin = init_2D_hist.FindBin(0,ylower)/init_2D_hist.GetNbinsX()
        # Range is inclusive, so don't count upper bin twice
        upper_bin = init_2D_hist.FindBin(0,yupper*(1-0.0001))/init_2D_hist.GetNbinsX()
        ybinned_hist = init_2D_hist.ProjectionX(hist_name, lower_bin, upper_bin, "e")
        ybinned_hist = ybinned_hist.Rebin(len(xbins)-1, hist_name+"_rebin", xbins)
        hists_half_rolled.append(ybinned_hist)

    if name is "":
        name = init_2D_hist.GetName().replace("2D", "unrolled")
    unrolled_hist = ROOT.TH1D(name, "Unrolled", nbins, 0, nbins)
    unrolled_hist.SetDirectory(init_2D_hist.GetDirectory())
    for i, hist in enumerate(hists_half_rolled):
        for j in range(1, hist.GetNbinsX()+1):
            entry = i*(hist.GetNbinsX()) + j
            content = hist.GetBinContent(j)
            error = hist.GetBinError(j)
            if j == hist.GetNbinsX():
                content += hist.GetBinContent(j+1)
                error += hist.GetBinError(j+1)
            unrolled_hist.SetBinContent(entry, content)
            unrolled_hist.SetBinError(entry, error)

    return unrolled_hist

def rebinHist(tmphist, histname, rebin):
    if rebin and type(tmphist) != ROOT.TH2D:
        if type(rebin) != int:
            return tmphist.Rebin(len(rebin)-1, histname, rebin)
        tmphist.Rebin(rebin)
    hist = tmphist.Clone(histname)
    ROOT.SetOwnership(hist, False)
    return hist

def make1DaQGCHists(orig_file, input2D_hists, plot_info, rebin=None):
    output_folders = []
    for name, data in plot_info.iteritems():
        entry = data["lheWeightEntry"]
        file_name = str(data["Members"][0])

        output_list = ROOT.TList()
        output_list.SetName(name)

        for init_2D_hist_name in input2D_hists:
            init_2D_hist = orig_file.Get("/".join([file_name, init_2D_hist_name]))
            # If a histogram with the same name exisits, ROOT will return
            # that instead of creating a new one. See:
            # https://root.cern.ch/root/html532/src/TH2.cxx.html#2253
            tmphist = init_2D_hist.ProjectionX("temphist", entry, entry, "e")
            histname = init_2D_hist_name.replace("lheWeights_", "")
            hist1D = rebinHist(tmphist, histname, rebin)
            tmphist.Delete()
            ROOT.SetOwnership(hist1D, False)
            output_list.Add(hist1D)
        output_folders.append(output_list)
    return output_folders

def removeZeros(hist):
    for i in range(hist.GetNbinsX()+2):
        if hist.GetBinContent(i) <= 0:
            if "Up" in hist.GetName():
                hist.SetBinContent(i, 0.0001)
            elif "Down" in hist.GetName():
                hist.SetBinContent(i, 0.00001)
            else: 
                hist.SetBinContent(i, 0.00005)

def getStatHists(hist, name, chan, signal):
    stat_hists = []
    variation_names = []
    for i in range(1, hist.GetNbinsX()+1):
        if "data" in name:
            continue
        if not (signal  == "wzjj_ewk" and "wzjj-vbfnlo" in name) \
                and not (signal  == "wzjj_vbfnlo" and "wzjj-ewk" in name):
            variation_names.append("%s_%s_statBin%i" % (name, chan, i))
        statUp_hist = hist.Clone(hist.GetName().replace(
            chan, "%s_%s_statBin%sUp_%s" % (name, chan, i, chan)))
        statDown_hist = hist.Clone(hist.GetName().replace(
            chan, "%s_%s_statBin%iDown_%s" % (name, chan, i, chan)))
        up = hist.GetBinContent(i)+hist.GetBinErrorUp(i) if \
                hist.GetBinContent(i) > 0 else hist.GetBinErrorUp(i)
        down = hist.GetBinContent(i)-hist.GetBinErrorLow(i)
        statUp_hist.SetBinContent(i, up) 
        statDown_hist.SetBinContent(i, down if down > 0 else 0.0001) 
        stat_hists.extend([statUp_hist, statDown_hist][:])
    for hist in stat_hists:
        removeZeros(hist)
    return (stat_hists, variation_names)

def getWeightHistProjection(init2D_hist, name, entry, rebin): 
    histname = init2D_hist.GetName().replace("lheWeights", name+"_weight%i" % entry)
    tmphist = init2D_hist.ProjectionX("temp", entry, entry, "e")
    hist = rebinHist(tmphist, histname, rebin)
    return hist

def getLHEWeightHists(init2D_hist, entries, name, variation_name, rebin=None):
    hists = []
    for i in entries:
        hist = getWeightHistProjection(init2D_hist, name, i, rebin)
        hists.append(hist)
    hist_name = init2D_hist.GetName()
    if "lheWeights" in hist_name:
        hist_name = hist_name.replace("lheWeights", "%s_%sUp" % (variation_name, name))
    else:
        hist_name = "_".join([hist_name, variation_name, name+"Up"])
    return hists, hist_name

def getMCPDFVariationHists(init2D_hist, entries, name, rebin=None, central=0):
    hists, hist_name = getLHEWeightHists(init2D_hist, entries, name, "pdfMC", rebin)
    if central == -1:
        upaction = lambda x: x[int(0.84*len(entries))] 
        downaction = lambda x: x[int(0.16*len(entries))] 
    else:
        upaction = lambda x : x[central]*(1+getPDFPercentVariation(x))
        downaction = lambda x: x[central]*(1-getPDFPercentVariation(x))

    return getVariationHists(hists, name, hist_name, 
            upaction, downaction, central
    )

# Calculate standard deviation per bin
def getSymmMCPDFVariationHists(init2D_hist, entries, name, rebin=None, central=0, pdfName=""):
    hists, hist_name = getLHEWeightHists(init2D_hist, entries, name, "pdf%sMC" % pdfName, rebin)
    upaction = lambda x: numpy.mean(x[1:]) + numpy.std(x[1:], ddof=1)
    downaction = lambda x: numpy.mean(x[1:]) - numpy.std(x[1:], ddof=1)

    return getVariationHists(hists, name, hist_name, 
            upaction, downaction, central
    )

def getAllSymmetricHessianVariationHists(init2D_hist, entries, name, rebin=None, central=0):
    hists, hist_name = getLHEWeightHists(init2D_hist, entries, name, "pdf", rebin)
    #centralIndex = central if central != -1 else int(len(entries)/2)
    variationSet = []
    for i, hist in enumerate(hists[0:central]+hists[central:]):
        upaction = lambda x: x[1] if x[1] > x[0] else (x[0]**2/x[1] if x[1] > 0 else 0)
        #downaction = lambda x: x[1] if x[1] < x[0] else float(x[0])/(x[1] if x[1] > 0 else 1)
        downaction = lambda x: x[1] if x[1] < x[0] else (x[0]**2/x[1] if x[1] > 0 else 0)
        variationSet.extend(getVariationHists([hists[central], hist], name, 
            hist_name.replace("pdf_%s" % name, "pdf%i" %i), upaction, downaction, central))
    return variationSet

def getHessianPDFVariationHists(init2D_hist, entries, name, rebin=None, central=0, pdfName=""):
    hists, hist_name = getLHEWeightHists(init2D_hist, entries, name, "pdf%sHes" % pdfName, rebin)
    sumsq = lambda x: math.sqrt(sum([0 if y < 0.01 else ((x[central] - y)**2) for y in x]))
    upaction = lambda x: x[central] + sumsq(x) 
    downaction = lambda x: x[central] - sumsq(x) 
    return getVariationHists(hists, name, hist_name, 
            upaction, downaction, central)

def getAssymHessianShift(vals, varType):
    pairs = zip(vals[1::2], vals[2::2])
    central = vals[0]
    diffs = [max(0, central - x[0], central - x[1]) for x in pairs] if varType == "down" else \
            [max(0, x[0] - central, x[1] - central) for x in pairs] 
    return math.sqrt(sum([x**2 for x in diffs]))

def getAssymHessianPDFVariationHists(init2D_hist, entries, name, rebin=None, central=0, pdfName=""):
    hists, hist_name = getLHEWeightHists(init2D_hist, entries, name, "pdf%sHes" % pdfName, rebin)
    #centralIndex = central if central != -1 else int(len(entries)/2)
    # From Eq. 3 https://arxiv.org/pdf/1101.0536.pdf
    #sumsqup = lambda x: math.sqrt(sum([max(0, i - x[central], j - x[central])**2 for i,j in zip(x[1::2], x[2::2])]))
    #sumsqdown = lambda x: math.sqrt(sum([max(0, x[central] - i, x[central] - j)**2 for i,j in zip(x[1::2], x[2::2])]))
    upaction = lambda x: x[central] + getAssymHessianShift(x, "up") 
    downaction = lambda x: x[central] - getAssymHessianShift(x, "down") 
    return getVariationHists(hists, name, hist_name, 
            upaction, downaction, central
    )

def getPDFPercentVariation(values):
    upvar = int(0.84*len(values))
    downvar = int(0.16*len(values))
    denom = values[upvar] + values[downvar]
    if denom == 0: 
        return 0
    return abs(values[upvar] - values[downvar])/denom

def getScaleHists(scale_hist2D, name, rebin=None, entries=[i for i in range(1,10)], central=0, exclude=[7,9]):
    hists, hist_name = getLHEWeightHists(scale_hist2D, entries, name, "QCDscale", rebin)
    return getVariationHists(hists, name, hist_name, lambda x: x[-1], lambda x: x[1], central)

# Pairs should correspond to muR, muF in (0.5, 1)
# nano ordering is [0] is mur=0.5 muf=0.5 ; [1] is mur=0.5 muf=1 ; [2] is mur=0.5 muf=2 ; [3] is mur=1 muf=0.5 ; 
# [4] is mur=1 muf=1 ; [5] is mur=1 muf=2 ; [6] is mur=2 muf=0.5 ; [7] is mur=2 muf=1 ; [8] is mur=2 muf=2 *
def getExpandedScaleHists(scale_hist2D, name, rebin=None, entries=[i for i in range(1,10)], central=-1, 
        pairs=[(1,7), (3,5), (0,8)]):
    hists, hist_name = getLHEWeightHists(scale_hist2D, entries, name, "QCDscale", rebin)
    variationSet = []
    for label, indices in zip(["muR", "muF", "muRmuF"], pairs):
        varhists = getVariationHists([hists[i] for i in indices], name,
            hist_name.replace("QCDscale","QCDscale_"+label), lambda x: x[-1], lambda x: x[1], -1)
        variationSet.extend(varhists)
    return variationSet

def getVariationHists(hists, process_name, histUp_name, up_action, down_action, central=0):
    histUp = hists[0].Clone(histUp_name)
    histDown = histUp.Clone(histUp_name.replace("Up", "Down"))
    
    histCentral = hists.pop(central).Clone() if central != -1 else None
    for i in range(0, histUp.GetNbinsX()+2):
        vals = []
        for hist in hists:
            vals.append(hist.GetBinContent(i))
        # This is a hack, but LHE weights were off by a factor of 2 for these samples
        if process_name in ["wlnu_jetbinned_nlo_cp5"]:
            vals = [x*2. for x in vals]
        vals.sort()
        vals.insert(0, histCentral.GetBinContent(i) if histCentral else 0)
        histUp.SetBinContent(i, up_action(vals))
        histDown.SetBinContent(i, down_action(vals))
        # For now, skip this check on aQGC for now, since they're screwed up
        if "aqgc" in process_name: continue
    logging.debug("For process %s, hist %s: Central, down, up: %s, %s, %s" % \
            (process_name, histUp_name, histCentral.Integral() if histCentral else 0, histDown.Integral(), histUp.Integral()))
    if histCentral and False: # Off for now, it can happen that groups have some hists with no weights which screws this up
        isValidVariation(process_name, histCentral, histUp, histDown)
    return [histUp, histDown] 

def isValidVariation(process_name, histCentral, histUp, histDown):
    for i in range(0, histCentral.GetNbinsX()+2):
        if histDown.GetBinContent(i) > histCentral.GetBinContent(i) and histCentral.GetBinContent(i) > 0.01:
            raise RuntimeError("Down variation >= central value for %s, hist %s"
                " This shouldn't be possible.\n"
                "up_hist: %0.4f\n" 
                "down_hist: %0.4f\n" 
                "central_hist: %0.4f\n" 
                "bin: %i\n" 
                % (process_name, histDown.GetName(), histUp.GetBinContent(i), histDown.GetBinContent(i), histCentral.GetBinContent(i), i)
            )
        if histUp.GetBinContent(i) < histCentral.GetBinContent(i) and histCentral.GetBinContent(i) > 0.01:
            raise RuntimeError("Up variation <= central value for %s, hist %s."
                " This shouldn't be possible.\n"
                "up_hist: %0.4f\n" 
                "down_hist: %0.4f\n" 
                "central_hist: %0.4f\n" 
                "bin: %i\n" 
                % (process_name, histUp.GetName(), histUp.GetBinContent(i), histDown.GetBinContent(i), histCentral.GetBinContent(i), i)
            )

def getTransformed3DScaleHists(scale_hist3D, transformation, transform_args, name):
    scale_hists = []
    for i in range(1,10):
        if i == 7 or i == 9: 
            continue
        scale_hist3D.GetZaxis().SetRange(i,i)
        # Order yx matters to have consistent axes!
        scale_hist2D = scale_hist3D.Project3D("yxe")
        scale_hist_name = scale_hist3D.GetName().replace("lheWeights", name+"_weight%i" % i)
        scale_hist2D.SetName(scale_hist_name)
        scale_hist1D = transformation(scale_hist2D, *transform_args)
        scale_hists.append(scale_hist1D)
    hist_name = scale_hist3D.GetName().replace("2D_lheWeights", "_".join(["unrolled", "QCDscale", name+"Up"]))
    return getVariationHists(scale_hists, name, hist_name, lambda x: x[-1], lambda x: x[1])

def getTransformed3DPDFHists(hist3D, transformation, transform_args, entries, name):
    hists = []
    for i in entries:
        hist3D.GetZaxis().SetRange(i,i)
        # Order yx matters to have consistent axes!
        hist2D = hist3D.Project3D("yxe")
        hist_name = hist3D.GetName().replace("lheWeights", name+"_weight%i" % i)
        hist2D.SetName(hist_name)
        hist1D = transformation(hist2D, *transform_args)
        hists.append(hist1D)
    #return hists
    hist_name = hist3D.GetName().replace("2D_lheWeights", "_".join(["unrolled", "pdf", name+"Up"]))
    return getVariationHists(hists, name, hist_name, 
            lambda x: x[0]*(1+getPDFPercentVariation(x)), 
            lambda x: x[0]*(1-getPDFPercentVariation(x))
    )

def addControlRegionToFitHist(control_hist, input_hist, base_name="unrolled"):
    hist = ROOT.TH1D("tmp", input_hist.GetTitle(), 
            input_hist.GetNbinsX()+1, 0, input_hist.GetNbinsX()+1)
    hist.SetName(input_hist.GetName().replace(base_name, base_name+"_wCR"))
    control_err = array.array('d', [0])
    control_yield = control_hist.IntegralAndError(0, control_hist.GetNbinsX()+1, control_err)
    hist.SetBinContent(1, control_yield) 
    hist.SetBinError(1, control_err[0])
    for i in range(1, hist.GetNbinsX()+1):
        hist.SetBinContent(i+1, input_hist.GetBinContent(i))
        hist.SetBinError(i+1, input_hist.GetBinError(i))
    ROOT.SetOwnership(hist, False)
    return hist

def addOverflow(hist):
    addOverflowAndUnderflow(hist, underflow=False, overflow=True)

def addOverflowAndUnderflow(hist, underflow=True, overflow=True):
    if not "TH1" in hist.ClassName():
        return
    if overflow:
        # Returns num bins + overflow + underflow
        num_bins = hist.GetNbinsX()
        add_overflow = hist.GetBinContent(num_bins) + hist.GetBinContent(num_bins + 1)
        hist.SetBinContent(num_bins, add_overflow)
    if underflow:
        add_underflow = hist.GetBinContent(0) + hist.GetBinContent(1)
        hist.SetBinContent(1, add_underflow)

def makeCompositeHists(hist_file, name, members, lumi, hists=[], underflow=False, overflow=True, rebin=None):
    composite = ROOT.TList()
    composite.SetName(name)
    for directory in [str(i) for i in members.keys()]:
        # For aQGC, the different plot groups should already be in their own files
        if "aqgc" in directory:
            directory = name
        
        # Support for negative cross section processes (for nonprompt)
        xseclookup = directory if directory in members.keys() or "-" in directory[0] else directory.split("__")[0]
        directory = directory.replace("-", "")
        if not hist_file.Get(directory):
            logging.warning("Skipping invalid filename %s" % directory)
            continue
        if hists == []:
            hists = [i.GetName() for i in hist_file.Get(directory).GetListOfKeys()]
        sumweights = 0
        isData = "data" in directory.lower() or "nonprompt" in directory.lower()
        if not isData:
            sumweights_hist = hist_file.Get("/".join([directory, "sumweights"]))
            if not sumweights_hist:
                logging.warning("Failed to find sumWeights for dataset %s" % directory)
            else:
                sumweights = sumweights_hist.Integral(1, sumweights_hist.GetNbinsX()+2)
                sumweights_hist.Delete()
        for histname in hists:
            if histname == "sumweights": continue
            tmphist = hist_file.Get("/".join([directory, histname]))
            if not tmphist: 
                raise RuntimeError("Failed to produce histogram %s" % "/".join([directory, histname]))
            toRebin = rebin and not "TH2" in tmphist.ClassName()
            hist = rebinHist(tmphist, histname, rebin)
            tmphist.Delete()
            if hist:
                sumhist = composite.FindObject(hist.GetName())
                xsec = members[xseclookup]
                if sumweights:
                    hist.Scale(xsec*1000*lumi/sumweights)
                elif not isData:
                    hist.Scale(1000*lumi*numpy.sign(xsec))
                    logging.warning("No sumWeights found for dataset %s, scaling only by lumi." % directory)
                addOverflowAndUnderflow(hist, underflow, overflow)
            else:
                raise RuntimeError("hist %s was not produced for "
                    "dataset %s!" % (histname, directory))
            if not sumhist:
                sumhist = hist.Clone()
                composite.Add(sumhist)
            else:
                sumhist.Add(hist)
            hist.Delete()
    return composite

def getTransformedHists(orig_file, folders, input_hists, transformation, transform_inputs):
    output_folders = []
    for folder in folders:
        output_list = ROOT.TList()
        output_list.SetName(folder)
        ROOT.SetOwnership(output_list, False)
        for input_hist_name in input_hists:
            orig_hist = orig_file.Get("/".join([folder, input_hist_name]))
            if not orig_hist:
                if "Fakes" not in input_hist_name and \
                    "Up" not in input_hist_name and "Down" not in input_hist_name:
                    logging.warning("Histogram %s not found for dataset %s. Skipping." % (input_hist_name, folder))
                continue
            new_hist = transformation(orig_hist, *transform_inputs)
            ROOT.SetOwnership(new_hist, False)
            output_list.Add(new_hist)
        output_folders.append(output_list)
    return output_folders 

def addaQGCTheoryHists(rtfile_name, plot_groups, base_hist_name):
    rtfile = ROOT.TFile(rtfile_name, "update")
    for name in plot_groups:
        if "__" not in name:
            continue
        aqgc_dir = rtfile.Get(name)
        aqgc_dir.cd()

        for chan in ["eee", "eem", "emm", "mmm"]:
            central_name = name.split("__")[0]
            varhist_name = "_".join([base_hist_name, "pdf_%sUp" % central_name, chan])
            hists = [varhist_name, varhist_name.replace("Up", "Down"), varhist_name.replace("pdf", "QCDscale"), 
                    varhist_name.replace("pdf", "QCDscale").replace("Up", "Down")]
            for hist_name in hists:
                base_hist = rtfile.Get("/".join([central_name, base_hist_name + "_" +chan]))
                aqgc_hist = rtfile.Get("/".join([name, base_hist_name + "_" +chan]))
                var_hist = rtfile.Get("/".join([central_name, hist_name]))
                aqgc_varhist = var_hist.Clone(hist_name.replace(central_name, name))
                for i in range(1, base_hist.GetNbinsX()+1):
                    if base_hist.GetBinContent(i) <= 0: 
                        continue
                    scale = aqgc_hist.GetBinContent(i)/base_hist.GetBinContent(i)
                    aqgc_varhist.SetBinContent(i, aqgc_varhist.GetBinContent(i)*scale)
                aqgc_varhist.Write()
                var_hist.Delete()
                base_hist.Delete()
                aqgc_hist.Delete()
        aqgc_dir.Delete()
    rtfile.Close()
