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
parser.add_argument("--mc2hes", action='store_true',
    help="Convert MC errors to hessian")
parser.add_argument("-c", "--central", type=str, default="wlnu_jetbinned_nlo_cp5",
    help="Sample to use as central value")
parser.add_argument("--files", type=lambda x: [i.strip() for i in x.split(",")], 
    default=[], help="Samples to add to output file")
parser.add_argument("-d", "--data", type=str, default="wlnu_nlo",
    help="Sample to use as dummy data")
parser.add_argument("-a", "--append", type=str, default="",
    help="Append to output folder name")
parser.add_argument("-b", "--fitvar", type=str, default="ptWmet",
    help="Variable to use in the fit")
parser.add_argument("-f", "--input_file", type=str, required=True,
    help="Input hist file")
parser.add_argument("-l", "--lumi", type=float, 
    default=35.9*0.7, help="lumi")
parser.add_argument("--noPdf", action='store_true', 
    help="don't add PDF uncertainties")
parser.add_argument("--noPtVSplit", action='store_true', 
    help="Don't split scale uncertainties by pt(V)")
parser.add_argument("-r", "--rebin", 
                    type=str, default=None, help="Rebin array: "
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
    "WGen/NanoAOD",
)

#plot_groups = ["wlnu_lo", "wlnu_lo_cp5", "wlnu_nlo", "wlnu_jetbinned_nlo", "wlnu_jetbinned_nlo_cp5", ]
plot_groups = args.files if args.files else ["wpmunu_minnlo_nnlopslike_photos", "wpmunu_nnlops_photos", "wpmunu_nnlops_nlow"]
plotGroupsMap = {name : config_factory.getPlotGroupMembers(name) for name in plot_groups}

xsecs  = ConfigureJobs.getListOfFilesWithXSec([f for files in plotGroupsMap.values() for f in files])

#channels = ["ep", "en", "mp", "mn"]
channels = ["mp", "mn"]
if args.rebin:
    if ":" in args.rebin:
        bins = UserInput.getRebin(args.rebin)
        args.rebin = array.array('d', bins)
    elif "," in args.rebin:
        args.rebin = array.array('d', [float(i.strip()) for i in args.rebin.split(",")])
    else:
        args.rebin = int(args.rebin)
    cardtool.setRebin(args.rebin)

cardtool.setFitVariable(args.fitvar)
if "unrolled" in args.fitvar:
    cardtool.setUnrolled([-2.5+0.5*i for i in range(0,11)], range(26, 56, 3))
cardtool.setProcesses(plotGroupsMap)
cardtool.setChannels(channels)
cardtool.setCrosSectionMap(xsecs)
cardtool.setVariations(["mWShift100MeV", "mWShift20MeV", "mWShift50MeV", "CMS_scale_m"])
folder_name = "_".join([args.fitvar,args.append]) if args.append != "" else args.fitvar
cardtool.setOutputFolder("/eos/user/k/kelong/CombineStudies/WGen/%s" % folder_name)

cardtool.setLumi(args.lumi)
cardtool.setInputFile(args.input_file)
cardtool.setOutputFile("WGenCombineInput.root")
#cardtool.setCombineChannels({"all" : channels, "e" : ["ep", "en"], "m" : ["mp", "mn"]})
#cardtool.setCombineChannels({"e" : ["ep", "en"], "m" : ["mp", "mn"]})
cardtool.setCombineChannels({"m" : ["mp"]})

ptbins = [0,3,5,7,9,12,15,20,27,40,100]
ptbinPairs = [(x,y) for x,y in zip(ptbins[:-1], ptbins[1:])]

for process in plot_groups:
    #Turn this back on when the theory uncertainties are added
    if "minnlo" in process:
        cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[7, 9], central=0)
        cardtool.setScaleVarGroups(process, [(3,6), (1,2), (4,8)])
        if not args.noPdf and "update" not in "minnlo":
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
        cardtool.setScaleVarGroups(process, [(3,6), (1,2), (4,8)])
        if not args.noPdf:
            # NNPDF3.1
            cardtool.addTheoryVar(process, 'pdf_hessian', range(885, 986), central=0, specName="NNPDF31")
    elif process not in ["nonprompt", "data"]:
        cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[3, 7], central=4)
        if not args.noPdf:
            cardtool.addTheoryVar(process, 'pdf_mc' if "cp5" not in process else "pdf_hessian", range(10,111), central=0)

    if not args.noPtVSplit:
        for pair in ptbinPairs:
            varName = 'ptV%ito%i' % pair
            varName = varName.replace("100", "Inf")
            cardtool.addScaleBasedVar(process, varName) 
    cardtool.addPerBinVariation(process, "CMS_eff_m", 0.01, False)

    cardtool.loadHistsForProcess(process, expandedTheory=True)
    cardtool.writeProcessHistsToOutput(process)

nuissance_map = {"e" : 105, "m" : 105 }
for chan in ["m"]:#["e", "m"]:
    cardtool.setTemplateFileName("Templates/CombineCards/VGen/WGen_template_{channel}.txt")
    logging.info("Writting cards for channel %s" % chan)
    cardtool.writeCards(chan, nuissance_map[chan], 
        extraArgs={"data_name" : args.data, 
            "w_sample" : args.central, "w_yield" : "yield:%s" % args.central, 
        }
    )

