# coding: utf-8
import ROOT
# Times bin size (0.5) and change units to pb
variables = {
    "ptZ" : { "histname" : "pT_lep1_lep2",
        "rebinScale" : 10*0.001,
    },
    "yZ" : { "histname" : "yZ",
        "rebinScale" : 0.1*0.001,
    },
}
pdf = "30"
outfile_name = "/eos/user/k/kelong/HistFiles/ZGen/NoSelection/ZToMuMu_MATRIX_NNPDF%s.root" % pdf

outfile = ROOT.TFile(outfile_name, "recreate")

for order in ["NLO", "NNLO", "NLO_NNLOPDF"]:
    file_path = "/eos/user/k/kelong/MatrixFiles/NNPDF%s" % pdf
    process = "DYm50_matrix" 
    if pdf == "30":
        process += "_nnpdf30" 
    if order[:3] == "NLO":
        process = process.replace("matrix", "matrix_nlo")
    if "NNLOPDF" in order:
        process += "_nnlopdf"
        file_path += "/NLO_NNLOPDF"
    outfile.mkdir(process)
    for key, value in variables.iteritems():
        name = "%s__%s_QCD" % (value["histname"], order.split("_")[0])
        print name
        rtfile = ROOT.TFile("%s/%s.root" % (file_path, name))
        hist = rtfile.Get("canvas").GetListOfPrimitives().FindObject(name)
        outfile.cd(process)
        for var in ["lhe", ""]:
            new_name = "_".join([key, "mm"] if var == "" else [key, var, "mm"])
            new_hist = hist.Clone(new_name)
            new_hist.Scale(value["rebinScale"]) 
            new_hist.Write()
