import ROOT

def float2double(hist):
    if hist.ClassName() == 'TH1D' or hist.ClassName() == 'TH2D':
        return hist
    elif hist.ClassName() == 'TH1F':
        new = ROOT.TH1D()
        hist.Copy(new)
    elif hist.ClassName() == 'TH2F':
        new = ROOT.TH2D()
        hist.Copy(new)
    else:
        raise Exception("Bad hist, dummy")
    return new

def invert2DHist(hist):
    new_hist = hist.Clone()
    ROOT.SetOwnership(new_hist, False)
    for x in range(hist.GetNbinsX()+1):
        for y in range(hist.GetNbinsY()+1):
            value = hist.GetBinContent(x, y)
            new_hist.SetBinContent(y, x, value)
    new_hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    new_hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    return new_hist

def setScaleFactorObj(infoVector, scaleFactorsObj, fScales):
    TotalSF = None
    for info in infoVector:
        # filename in 0, histname in 0
        f = ROOT.TFile(info[0])
        hist = float2double(f.Get(info[1]))
        if TotalSF == None:
            TotalSF = hist.Clone()
        else:
            TotalSF.Multiply(hist)
    scaleFactorsObj.Set2DHist(TotalSF)
    fScales.cd()
    scaleFactorsObj.Write()
    

