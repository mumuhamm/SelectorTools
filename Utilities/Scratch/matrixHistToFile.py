# coding: utf-8
import ROOT
chan = "mn"
resum = True
fiducial = True
variables = {
    "ptZ" : { "histname" : "pT_lep1_lep2",
        "rebinScale" : 0.001,
    },
    "yZ" : { "histname" : "yZ",
        "rebinScale" : 0.1*0.001,
    },
    "ptl1" : { "histname" : "pT_lep1",
        "rebinScale" : 2*0.001,
    },
    "ptl2" : { "histname" : "pT_lep2",
        "rebinScale" : 2*0.001,
    },
}
pdf = "31"
outfile_name = "/eos/user/k/kelong/HistFiles/ZGen/NoSelection/ZToMuMu_MATRIX_EWParamMatch_NNPDF%s.root" % pdf

if chan != "mm":
    outfile_name = "/eos/user/k/kelong/HistFiles/WGen/NoSelection/WpToMuNu_MATRIX_EWParamMatch_NNPDF%s.root" % pdf
    #outfile_name = "/eos/user/k/kelong/HistFiles/WGen/NoSelection/WpToMuNu_MATRIX_DefaultParams_NNPDF%s.root" % pdf
    #outfile_name = "/eos/user/k/kelong/HistFiles/WGen/NoSelection/WpToMuNu_MATRIX_Marius.root"
    if chan == "mn":
        outfile_name = outfile_name.replace("Wp", "Wm")
    variables = {
        "ptW" : { "histname" : "pT_W",
            "rebinScale" : 1*0.001,
        },
        "yW" : { "histname" : "yW",
            "rebinScale" : 0.05*0.001,
        },
        "ptl" : { "histname" : "pT_lep",
            "rebinScale" : 0.001,
        },
        "ptnu" : { "histname" : "pT_nu",
            "rebinScale" : 0.001,
        },
    }

if resum:
    if chan == "mm":
        variables = {
            "ptZ" : { "histname" : "RadISH_observable_binning2",
                "rebinScale" : 0.001,
            },
        }
        outfile_name = "/eos/user/k/kelong/HistFiles/ZGen/NoSelection/ZToMuMu_MATRIX_RadISH_MatchEWParams_NNPDF%s.root" % pdf
    else:
        variables = {
            "ptW" : { "histname" : "RadISH_observable_binning2",
                "rebinScale" : 0.001,
            },
        }
        outfile_name = "/eos/user/k/kelong/HistFiles/WGen/NoSelection/WmToMuNu_MATRIX_RadISH_MatchEWParams_NNPDF%s.root" % pdf


if fiducial:
    outfile_name = outfile_name.replace("/NoSelection", "")
outfile = ROOT.TFile(outfile_name, "recreate")

#for order in ["NLO", "NNLO", "NLO_NNLOPDF"]:
for order in ["NNLO", ]:
    file_path = "/eos/user/k/kelong/MatrixFiles/NNPDF%s_MatchEWParams" % pdf
    if fiducial:
        file_path += "/Fiducial"

    process = "DYm50_matrix" 
    if chan == "mp":
        process = "wpmunu_matrix"
    elif chan == "mn":
        process = "wmmunu_matrix"
        file_path = "/".join([file_path, "Wm"])

    elif "Default" in outfile_name:
        file_path = "/eos/user/k/kelong/MatrixFiles/NNPDF%s" % pdf
        process = process+"__default" 
    elif "Marius" in outfile_name and chan != "mm":
        file_path = "/eos/user/k/kelong/MatrixFiles/Marius" 
        process = process+"__marius" 
        variables["yW"]["histname"] = "yV"
        variables["yW"]["rebinScale"] = 0.1*0.001
        variables["ptW"]["histname"] = "ptV"

    if pdf == "30":
        process += "_nnpdf30" 
    if order[:3] == "NLO":
        process = process.replace("matrix", "matrix_nlo")
    if "NNLOPDF" in order:
        process += "_nnlopdf"
        file_path += "/NLO_NNLOPDF"
    outfile.mkdir(process)
    for key, value in variables.iteritems():
        proc = process
        if resum and "pt" in key:
            proc = process + "__radish"
            outfile.mkdir(proc)
            file_path = "/eos/user/k/kelong/MatrixFiles/Radish"
            if chan != "mm" and fiducial:
                file_path = file_path + "/Wm_Fiducial"
            order += "+N3LL"
        print file_path
        name = "%s__%s_QCD" % (value["histname"], order.split("_")[0])
        if resum:
            name = name.replace("_QCD", "")
        rtfile = ROOT.TFile("%s/%s.root" % (file_path, name))
        hist = rtfile.Get("canvas").GetListOfPrimitives().FindObject(name)
        uphist = rtfile.Get("canvas").GetListOfPrimitives().FindObject(name+"__scaleUp")
        downhist = rtfile.Get("canvas").GetListOfPrimitives().FindObject(name+"__scaleDown")
        outfile.cd(proc)
        basevars = ["QCDscale_"+proc+"Up", "QCDscale_"+proc+"Down",]
        variations = basevars + ["", "lhe"] + ["lhe_"+v for v in basevars]
        
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
            empty_hist = tmp.Clone(new_name.replace(chan, chan.replace("m", "e")))
            empty_hist.Scale(0)
            empty_hist.Write()
print outfile
