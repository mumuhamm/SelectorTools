from python import ConfigureJobs,CombineCardTools,UserInput
import sys
import ROOT
import logging
import array
import argparse

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("--debug", action='store_true',
    help="Print debug info")
parser.add_argument("-a", "--append", type=str, default="",
    help="Append to output folder name")
parser.add_argument("-f", "--input_file", type=str, required=True,
    help="Input hist file")
parser.add_argument("-b", "--fitvar", 
    type=str, default="ptZ",
    help="Variable to use in the fit")
parser.add_argument("-c", "--central", type=str, default="wlnu_jetbinned_nlo_cp5",
    help="Sample to use as central value")
parser.add_argument("-d", "--data", type=str, default="wlnu_nlo",
    help="Sample to use as dummy data")
parser.add_argument("--noPdf", action='store_true', 
    help="don't add PDF uncertainties")
parser.add_argument("--files", type=lambda x: [i.strip() for i in x.split(",")], 
    default=[], help="Samples to add to output file")
parser.add_argument("-r", "--rebin", 
                    type=str, default=None, help="Rebin array: "
                    "values (bin edges) separated by commas.")
args = parser.parse_args()

channels = ["mp", "mn"]
logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO))

cardtool = CombineCardTools.CombineCardTools()

if args.rebin:
    if ":" in args.rebin:
        args.rebin = array.array('d', UserInput.getRebin(args.rebin))
    elif "," in args.rebin:
        args.rebin = array.array('d', [float(i.strip()) for i in args.rebin.split(",")])
    else:
        args.rebin = int(args.rebin)
    cardtool.setRebin(args.rebin)

manager_path = ConfigureJobs.getManagerPath() 
sys.path.append("/".join([manager_path, "AnalysisDatasetManager",
    "Utilities/python"]))

from ConfigHistFactory import ConfigHistFactory
config_factory = ConfigHistFactory(
    "%s/AnalysisDatasetManager" % manager_path,
    "ZGen/NanoAOD",
)

plot_groups = args.files if args.files else ["dy_minnlo_photos", "dy_nlo", ]
#plot_groups = ["dy_htbinned_cp5", "dy_lo_cp5", "dy_nlo_cp5", "dy_lo", "dy_nlo", "dy_nlo_jetbinned", "dy_nlo_jetbinned_cp5"]
plotGroupsMap = {name : config_factory.getPlotGroupMembers(name) for name in plot_groups}

xsecs  = ConfigureJobs.getListOfFilesWithXSec([f for files in plotGroupsMap.values() for f in files])

#channels = ["ee", "mm"]
channels = ["mm"]
rebin = array.array('d', args.rebin) if args.rebin else None
cardtool.setFitVariable(args.fitvar)
if rebin:
    cardtool.setRebin(rebin)
cardtool.setProcesses(plotGroupsMap)
cardtool.setChannels(channels)
cardtool.setCrosSectionMap(xsecs)
cardtool.setVariations([])
folder_name = "_".join([args.fitvar,args.append]) if args.append != "" else args.fitvar
cardtool.setOutputFolder("/eos/user/k/kelong/CombineStudies/ZGen/%s" % folder_name)

cardtool.setLumi(35.9)
cardtool.setInputFile(args.input_file)
cardtool.setOutputFile("ZGenCombineInput.root")

ptbins = [0,3,5,7,9,12,15,20,27,40,100]
ptbinPairs = [(x,y) for x,y in zip(ptbins[:-1], ptbins[1:])]

for process in plot_groups:
    #Turn this back on when the theory uncertainties are added
    if "update_ref" in process or "lowcutoff" in process:
        cardtool.addTheoryVar(process, 'scale', range(10, 19), exclude=[15, 17], central=0)
    elif "minnlo" in process:
        cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[6, 8], central=0)
        if not args.noPdf:
            # NNPDF3.1
            cardtool.addTheoryVar(process, 'pdf_hessian', range(10, 111), central=0, specName="NNPDF31")
            # NNPDF31_nnlo_as_0118_CMSW1_hessian_100; LHAPDFID = 325700
            cardtool.addTheoryVar(process, 'pdf_hessian', range(121, 222), central=0, specName="CMSW1")
            # NNPDF31_nnlo_as_0118_CMSW2_hessian_100; LHAPDFID = 325900
            cardtool.addTheoryVar(process, 'pdf_hessian', range(222, 323), central=0, specName="CMSW2")
            # NNPDF31_nnlo_as_0118_CMSW3_hessian_100; LHAPDFID = 326100
            cardtool.addTheoryVar(process, 'pdf_hessian', range(323, 424), central=0, specName="CMSW3")
            # NNPDF31_nnlo_as_0118_CMSW3_hessian_100; LHAPDFID = 326300
            cardtool.addTheoryVar(process, 'pdf_hessian', range(424, 525), central=0, specName="CMSW4")
            # CT14
            cardtool.addTheoryVar(process, 'pdf_assymhessian', range(525, 582), central=0, specName="CT14")
            # MMHT
            cardtool.addTheoryVar(process, 'pdf_assymhessian', range(584, 635), central=0, specName="MMHT")
            # HERA20_EIG
            cardtool.addTheoryVar(process, 'pdf_assymhessian', range(668, 711), central=0, specName="HERA2")
    elif "nnlops" in process:
        cardtool.addTheoryVar(process, 'scale', range(10, 19), exclude=[15, 17], central=0)
    elif process not in ["nonprompt", "data"]: #and False
        cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[6, 7], central=4)
        pdf_entries = [4] + (range(10, 40) if "cp5" in process else range(10, 110))
        #cardtool.addTheoryVar(process, 'pdf_mc' if "cp5" in process else "pdf_hessian", pdf_entries, central=0)
        cardtool.addTheoryVar(process, 'pdf_hessian', pdf_entries, central=0)

    for pair in ptbinPairs:
        varName = 'ptV%ito%i' % pair
        varName = varName.replace("100", "Inf")
        cardtool.addScaleBasedVar(process, varName) 

    cardtool.loadHistsForProcess(process)
    cardtool.writeProcessHistsToOutput(process)

nuissance_map = {"ee" : 5, "mm" : 5 }
for chan in channels: #+ ["all"]:
    cardtool.setTemplateFileName("Templates/CombineCards/VGen/ZGen_template_{channel}.txt")
    logging.info("Writting cards for channel %s" % chan)
    cardtool.writeCards(chan, nuissance_map[chan], 
        extraArgs={"data_name" : args.data, 
            "dy_sample" : args.central, "dy_yield" : "yield:%s" % args.central,
        }
    )
