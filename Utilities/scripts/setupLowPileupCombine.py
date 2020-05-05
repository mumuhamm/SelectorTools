from python import CombineCardTools
from python import ConfigureJobs
import sys
import ROOT
import logging
import array
import argparse

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("--debug", action='store_true',
    help="Print debug info")
parser.add_argument("-u", '--unbinnedSignal', action='store_true',
    help="Just use one signalStrength")
parser.add_argument("-a", "--append", type=str, default="",
    help="Append to output folder name")
parser.add_argument("-f", "--input_file", type=str, required=True,
    help="Input hist file")
parser.add_argument("-b", "--fitvars", 
    type=lambda x: x.split(","), default=["ptW"],
    help="Variable to use in the fit")
parser.add_argument("-r", "--rebin", 
                    type=lambda x : [] if "," not in x else [float(i.strip()) for i in x.split(',')],
                    default=[i for i in range(0, 105, 5)], help="Rebin array: "
                    "values (bin edges) separated by commas.")
args = parser.parse_args()

logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO))

cardtool = CombineCardTools.CombineCardTools()

manager_path = ConfigureJobs.getManagerPath() 
sys.path.append("/".join([manager_path, "AnalysisDatasetManager",
    "Utilities/python"]))

from ConfigHistFactory import ConfigHistFactory
config_factory = ConfigHistFactory(
    "%s/AnalysisDatasetManager" % manager_path,
    "LowPileupW/NanoAOD",
)

plot_groups = ["data", "nonprompt_m", "vv", "dy", "wtv", "top", "wlnu_jetbinned_nlo", ]

#ptbins = [0.0, 13.0, 26.0, 38.0, 50.0, 62.0, 75.0, 100.0]
if not args.unbinnedSignal:
    ptbins = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100]
    regions = ["pt%ito%i" % (ptbins[i], ptbins[i+1]) for i in range(len(ptbins)-1)]
    plot_groups += ["wmv_jetbinned_nlo_%s" % region for region in regions]
    plot_groups.append("wmv_jetbinned_nlo_pt100")

plotGroupsMap = {name : config_factory.getPlotGroupMembers(name) for name in plot_groups}
print plotGroupsMap

#xsecs  = ConfigureJobs.getListOfFilesWithXSec([f for files in plotGroupsMap.values() for f in files])
files = ConfigureJobs.getListOfFiles([f for files in plotGroupsMap.values() for f in files], selection="Wselection")
xsecs = { f : 1. for f in files }

channels = ["mp", "mn"]
if args.rebin:
    rebin = array.array('d', args.rebin)
    cardtool.setRebin(rebin)
cardtool.setProcesses(plotGroupsMap)
cardtool.setChannels(channels)
cardtool.setCrosSectionMap(xsecs)
variations = ["CMS_eff_MCsubt_m", "CMS_modeling_fsr", "CMS_eff_stat_m", 
        "CMS_eff_background_m", "CMS_eff_tagPt_m",
        "CMS_recoilCorrection_EtaShape",
        "CMS_recoilCorrection_RUShape", 
        "CMS_recoilCorrection_KeysShape",
        "CMS_recoilCorrection_StatBin0",    
        "CMS_recoilCorrection_StatBin1",    
        "CMS_recoilCorrection_StatBin2",        
        "CMS_recoilCorrection_StatBin3",        
        "CMS_recoilCorrection_StatBin4",        
        "CMS_recoilCorrection_StatBin5",        
        "CMS_recoilCorrection_StatBin6",        
        "CMS_recoilCorrection_StatBin7",        
        "CMS_recoilCorrection_StatBin8",        
        "CMS_recoilCorrection_StatBin9",        
]
cardtool.setVariations(variations, exclude=["nonprompt_e", "nonprompt_m", "data"])
folder_name = "_".join(args.fitvars + ([args.append] if args.append != "" else []))
cardtool.setOutputFolder("/eos/user/k/kelong/CombineStudies/LowPileup/%s" % folder_name)

cardtool.setFitVariableAppend("nonprompt_m", "Fakes")
cardtool.setFitVariableAppend("nonprompt_e", "Fakes")

cardtool.setLumi(0.199)
cardtool.setInputFile(args.input_file)
for fitvar in args.fitvars:
    cardtool.setFitVariable(fitvar)
    cardtool.setOutputFile("WLowPileup_%s_CombineInput.root" % fitvar)
    cardtool.setCombineChannels({"all" : channels, "m" : ["mp", "mn"]})
    for process in plot_groups:
        #Turn this back on when the theory uncertainties are added
        if process not in ["nonprompt_m", "nonprompt_e", "data"] and False:
            cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[6, 7], central=4)
            if "cp5" not in process and args.mc2hes:
                # Exclude alpha_s variations
                pdf_entries = [5] + range(112, 172) 
                cardtool.addTheoryVar(process, "pdf_hessian", pdf_entries, central=0)
            else:
                pdf_entries = [5] + (range(10, 110) if "cp5" not in process else range(10, 40))
                cardtool.addTheoryVar(process, 'pdf_mc' if "cp5" not in process else "pdf_hessian", pdf_entries, central=0)
        cardtool.loadHistsForProcess(process)
        cardtool.writeProcessHistsToOutput(process)
    cardtool.writeMetaInfo()

    nuissance_map = {"e" : len(variations)+6, "m" : len(variations)+6}
    template_name = "W" if args.unbinnedSignal else "WPtBinned"
    for chan in ["m"]:
        cardname = "%s_template_{channel}.txt" % template_name
        cardtool.setTemplateFileName("Templates/CombineCards/LowPileup/"+cardname)
        logging.info("Writting cards for channel %s" % chan)
        cardtool.writeCards(chan, nuissance_map[chan], outlabel=fitvar)
