# coding: utf-8
import ROOT
w = True
variables = {
    "ptZ" : { "histname" : "pT_lep1_lep2",
        "rebinScale" : 10*0.001,
    },
    "yZ" : { "histname" : "yZ",
        "rebinScale" : 0.1*0.001,
    },
    #"ptl1" : { "histname" : "pT_lep1",
    #    "rebinScale" : 0.1*0.001,
    #},
}
pdf = "31"
outfile_name = "/eos/user/k/kelong/HistFiles/ZGen/NoSelection/ZToMuMu_MATRIX_EWParamMatch_NNPDF%s.root" % pdf
chan="mm"

if w:
    outfile_name = "/eos/user/k/kelong/HistFiles/WGen/NoSelection/WpToMuNu_MATRIX_EWParamMatch_NNPDF%s.root" % pdf
    variables = {
        "ptW" : { "histname" : "pT_W",
            "rebinScale" : 10*0.001,
        },
        "yW" : { "histname" : "yW",
            "rebinScale" : 0.5*0.001,
        },
        #"ptl1" : { "histname" : "pT_lep1",
        #    "rebinScale" : 0.1*0.001,
        #},
    }
    chan="mp"

outfile = ROOT.TFile(outfile_name, "recreate")

#for order in ["NLO", "NNLO", "NLO_NNLOPDF"]:
for order in ["NNLO", ]:
    file_path = "/eos/user/k/kelong/MatrixFiles/NNPDF%s_MatchEWParams" % pdf
    process = "DYm50_matrix" if not w else "wpmunu_matrix"
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
        uphist = rtfile.Get("canvas").GetListOfPrimitives().FindObject(name+"__scaleUp")
        downhist = rtfile.Get("canvas").GetListOfPrimitives().FindObject(name+"__scaleDown")
        outfile.cd(process)
        basevars = ["QCDscale_"+process+"Up", "QCDscale_"+process+"Down",]
        variations = basevars + ["lhe"] + ["lhe_"+v for v in basevars]
        
        for var in variations:
            new_name = "_".join([key, chan] if var == "" else [key, var, chan])
            tmp = hist
            if "Up" in var:
                tmp = uphist
            elif "Down" in var:
                tmp = downhist
            new_hist = tmp.Clone(new_name)
            new_hist.Scale(value["rebinScale"]) 
            new_hist.Write()
