from python import HistTools, OutputTools, ConfigureJobs
weight id='1950'> lhapdf=263000 </weight>
import ROOT
import argparse
import json
import array

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", type=str,
    default="test.root", help="Output file name")
parser.add_argument("--rebin", type=lambda x: [float(i) for i in x.split(",")], 
                    default=ConfigureJobs.getBinning(),
                    help="Rebin values, comma separated list")
args = parser.parse_args()

rebin = array.array('d', args.rebin)
saveToFile = True
input_file_name = args.input_file
input_file = ROOT.TFile(input_file_name, "update" if saveToFile else "read")

chans = ["eee", "eem", "emm", "mmm"]
plot_group_info = {}
for file_name in ["WZxsec2016_aQGC-" + i for i in ["FT", "FM", "FS"]]:
    with open(ConfigureJobs.getManagerPath() + 
            "/AnalysisDatasetManager/PlotGroups/%s.json" % file_name) as plot_file:
        plot_group_info.update(json.load(plot_file))

jevariations = [i for x in ["CMS_scale_j", "CMS_res_j",] for i in [x+"Up", x+"Down"]]
variations = [i for x in ["CMS_scale_j", "CMS_res_j", \
    "CMS_eff_m", "CMS_scale_m", "CMS_eff_e", "CMS_scale_e", "CMS_pileup", "CMS_scale_unclEnergy"] for i in [x+"Up", x+"Down"]]
aqgc_hists = HistTools.make1DaQGCHists(input_file, 
        ["mjj_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
            ["Mass_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
            ["m_l1l3_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
            ["m_l2l3_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
        #    ["Pt_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
        #    ["llPt_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
        #    ["l2Pt_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
        #    ["l3Pt_lheWeights_%s" % c for c in ConfigureJobs.getChannels()] + \
            ["MTWZ_lheWeights_%s_%s" % (var, c) for c in ConfigureJobs.getChannels()
                for var in variations] + \
            ["MTWZ_lheWeights_%s" % c for c in ConfigureJobs.getChannels()],
        plot_group_info, rebin,
)

if saveToFile:
    for output in aqgc_hists:
        OutputTools.writeOutputListItem(output, input_file)

addControlRegion=True

if addControlRegion:
    for folder_name in input_file.GetListOfKeys():
        if folder_name.GetName() == "NonpromptMC":
            continue
        hist_names = []
        if "nonprompt" in folder_name.GetName().lower() or folder_name.GetName() == "DataEWKCorrected":
            hist_names = ["MTWZ_%s_Fakes_%s" % (var, chan) for var in jevariations for chan in chans]
            hist_names.extend(["MTWZ_Fakes_%s" % chan for chan in chans])
        else: 
            hist_names = ["MTWZ_%s_%s" % (var, chan) for var in variations for chan in chans]
            hist_names.extend(["MTWZ_%s" % chan for chan in chans])
        folder = input_file.Get(folder_name.GetName())
        folder.cd()
        for hist_name in hist_names:
            tmphist = folder.Get(hist_name)
            h = tmphist.Clone() if not rebin else tmphist.Rebin(len(rebin)-1, "", rebin)
            append = h.GetName().split("_")
            append = append[append.index("MTWZ")+1:]
            hist_name = "_".join(["backgroundControlYield"] + append)
            base_folder = folder.GetName().split("__")[0]
            control_hist = input_file.Get(base_folder + "/" +hist_name)
            if not control_hist:
                raise RuntimeError("failed to find hist %s for folder %s" % (hist_name, folder))
            hist = HistTools.addControlRegionToFitHist(control_hist, h, "MTWZ")
            hist.Write()
            h.Delete()
        input_file.cd()

