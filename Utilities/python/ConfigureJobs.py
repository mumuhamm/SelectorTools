import datetime
import UserInput
import fnmatch
import glob
import os
import json
import array
import string
import configparser

def get2DBinning(xvar="mjj", yvar="etajj", analysis='WZ'):
    #return (array.array('d', [500, 1000,1500, 2000, 2500]),
    # [0, 150, 300, 450] # for MT(WZ)
#    return (array.array('d', [500, 1000, 1350, 1750, 2000, 2500]),
    xbinning = []
    ybinning = []
    if xvar == "mjj":
        xbinning = array.array('d', [500, 1000,1500, 2000, 2500])
        #xbinning = array.array('d', [500, 1000, 1350, 1750, 2500])

    if yvar == 'etajj':    
        ybinning = [2.5, 4, 5, 20]
    #if yvar == 'etajj':    
    #    ybinning = [2.5, 4, 5.5, 20]
    elif yvar == 'dRjj':
        ybinning = [0, 5, 6, 20]
    return (xbinning, ybinning)

def getChannels(analysis='ZZ'):
    if analysis == 'ZZ':
        return ["eeee", "eemm", "mmee", "mmmm"]
def getManagerPath():
    config = configparser.ConfigParser()
    config.read_file(open("Templates/config.%s" % os.environ["USER"]))
    if "dataset_manager_path" not in config['Setup']:
        raise ValueError("dataset_manager_path not specified in config file Template/config.%s" 
                            % os.environ["USER"])
    return config['Setup']['dataset_manager_path'] + "/"
def getCombinePath():
    config = configparser.ConfigParser()
    config.read_file(open("Templates/config.%s" % os.environ["USER"]))
    if "combine_path" not in config['Setup']:
        raise ValueError("dataset_manager_path not specified in config file Template/config.%s" 
                            % os.environ["USER"])
    return config['Setup']['combine_path'] + "/"
def getHarvesterPath():
    config = configparser.ConfigParser()
    config.read_file(open("Templates/config.%s" % os.environ["USER"]))
    if "harvester_path" not in config['Setup']:
        raise ValueError("dataset_manager_path not specified in config file Template/config.%s" 
                            % os.environ["USER"])
    return config['Setup']['harvester_path'] + "/"
def getListOfEWKFilenames():
    return [
        "wz3lnu-mg5amcnlo",
    # Use jet binned WZ samples for subtraction by default
        #"wz3lnu-mgmlm-0j",
        #"wz3lnu-mgmlm-1j",
        #"wz3lnu-mgmlm-2j",
        #"wz3lnu-mgmlm-3j",
        #"wlljj-ewk",
        "zz4l-powheg",
        #"zz4ljj-ewk",
        #"zz2l2vjj-ewk",
        #"tzq",
        #"ttz",
        #"ttw",
        #"zzz",
        #"wwz",
        #"www",
        #"ww",
        #"zg",
        "ggZZ4e",
        "ggZZ4m",
        "ggZZ4t",
        "ggZZ2e2mu",
        "ggZZ2e2tau",
        "ggZZ2mu2tau",
    ]

def getListOfVVV():
    return [
       "ZZGJetsTo4L2Nu",
       "ZZZJetsTo4L2Nu",
       "WZZJetsTo4L2Nu",
       "WWZJetsTo4L2Nu", 
    ]
def getListOfEWK():
    return [
        "zz4l-powheg",
        "ggZZ4e",
        "ggZZ4m",
        "ggZZ4t",
        "ggZZ2e2mu",
        "ggZZ2e2tau",
        "ggZZ2mu2tau",
    ]
def getListOfSignalFilenames():
    return [
        "zz4l-powheg",
    ]
def getListOfDYFilenames():
    return[
        "DYJetsToLL_M10to50",
        "DYJetsToLLM-50",
        ]
def getListOfNonpromptFilenames():
    return ["tt-jets",
        #"st-schan",
        #"st-tchan-t",
        #"st-tchan-tbar",
        #"st-tw",
        #"st-tbarw",
        ##"DYm50",
        #"DYm50-1j",
        #"DYm50-2j",
        #"DYm50-3j",
        #"DYm50-4j",
    ]
def getListOfHZZFilenames():
    return ["ggHZZ",
            "vbfHZZ",
            "ttH_HToZZ_4L",
            "WminusHToZZ",
            "WplusHToZZ",
            "ZHToZZ_4L"
    ]
def getListOfggZZFilenames():
    return ["ggZZ4e",
        "ggZZ4m",
        "ggZZ4t",
        "ggZZ2e2mu",
        "ggZZ2e2tau",
        "ggZZ2mu2tau"
        ]
def getJobName(sample_name, analysis, selection, version):
    date = '{:%Y-%m-%d}'.format(datetime.date.today())
    selection = selection.replace(";",",")
    selections = selection.split(",")
    selection_name = "To".join([selections[0],selections[-1]]) \
        if len(selections) > 1 else selections[0]
    return '-'.join([date, sample_name, analysis, selection_name, 
        ("v%s" % version) if version.isdigit() else version])
def getNumberAndSizeOfLocalFiles(path_to_files):
    file_list = glob.glob(path_to_files)
    return (len(file_list), sum([os.path.getsize(f)/1000000 for f in file_list]))
def getNumberAndSizeOfHDFSFiles(file_path):
    file_info = []
    for line in out.splitlines():
        split = line.split()
        if len(split) != 9:
            continue
        file_info.append(float(split[4].rstrip("mkg")))
    return (0,0) if len(file_info) == 0 else (len(file_info), sum(file_info))
def getListOfHDFSFiles(file_path):
    try:
        out = subprocess.check_output(["hdfs", "dfs", "-ls", file_path.replace("/hdfs", "")])
    except subprocess.CalledProcessError as error:
        logging.warning(error)
        return []
    files = []
    for line in out.splitlines():
        split = line.split(" ", 1)
        if len(split) != 2:
            continue
        elif "root" in split[1]:
            files.append("/"+split[1])
    return files

def getMCInfo(filelist, selection, manager_path=""):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/AnalysisDatasetManager/FileInfo" % manager_path
    data_info = UserInput.readAllJson("/".join([data_path, "%s.json" % "data/*"]))
    mc_info = UserInput.readAllJson("/".join([data_path, "%s.json" % "montecarlo/*"]))
    
    return mc_info
def getListOfFiles(filelist, selection, manager_path=""):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/AnalysisDatasetManager/FileInfo" % manager_path
    data_info = UserInput.readAllJson("/".join([data_path, "%s.json" % "data/*"]))
    mc_info = UserInput.readAllJson("/".join([data_path, "%s.json" % "montecarlo/*"]))
    valid_names = data_info.keys() + mc_info.keys()
    names = []
    #print "name in getListOfFiles function: ",name
    for name in filelist:
        if "ZZ4l2018" in name:
            dataset_file = manager_path + \
                "AnalysisDatasetManager/FileInfo/ZZ4l2018/%s.json" % selection
            allnames = json.load(open(dataset_file)).keys()
            print allnames
            if "nodata" in name:
                nodata = [x for x in allnames if "data" not in x]
                names += nodata
            elif "data" in name:
                names += [x for x in allnames if "data" in x]
            else:
                names += allnames
        elif "ZZ4l2019" in name:
            dataset_file = manager_path + \
                "AnalysisDatasetManager/FileInfo/ZZ4l2019/%s.json" % selection
            allnames = json.load(open(dataset_file)).keys()
            print allnames
            if "nodata" in name:
                nodata = [x for x in allnames if "data" not in x]
                names += nodata
            elif "data" in name:
                names += [x for x in allnames if "data" in x]
            else:
                names += allnames
        elif "zz4l-powheg" in name:
            names += [name]
        elif "*" in name:
            names += fnmatch.filter(valid_names, name)
        else:
            if name.split("__")[0] not in valid_names:
                print "%s is not a valid name" % name
                continue
            names += [name]
    return [str(i) for i in names]
def fillTemplatedFile(template_file_name, out_file_name, template_dict):
    with open(template_file_name, "r") as templateFile:
        source = string.Template(templateFile.read())
        result = source.substitute(template_dict)
    with open(out_file_name, "w") as outFile:
        outFile.write(result)
def getListOfFilesWithXSec(filelist, manager_path=""):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/AnalysisDatasetManager/FileInfo" % manager_path
    files = getListOfFiles(filelist, "ntuples", manager_path)
    #files = getListOfFiles(filelist, "2018Data", manager_path)
    mc_info = UserInput.readAllJson("/".join([data_path, "%s.json" % "montecarlo/*"]))
    info = {}
    for file_name in files:
        if "data" in file_name:
            info.update({file_name : 1})
        else:
            file_info = mc_info[file_name.split("__")[0]]
            kfac = file_info["kfactor"] if "kfactor" in file_info.keys() else 1
            info.update({file_name : file_info["cross_section"]*kfac})
    return info
def getPreviousStep(selection, analysis):
    selection_map = {}
    if analysis == "ZZ4l2018":
        selection_map = { "ntuples" : "ntuples",
                "preselection" : "ntuples",
                "4lCRBase" : "ntuples"
        }
    elif analysis == "ZplusL2018":
        selection_map = { "ntuples": "ntuples",
                "ZplusLBase" : "ntuples"
        }
    elif analysis == "ZZ4l2019":
        selection_map = { "ntuples" : "ntuples",
                "preselection" : "ntuples",
                "4lCRBase" : "ntuples"
        }
    elif analysis == "ZplusL2019":
        selection_map = { "ntuples": "ntuples",
                "ZplusLBase" : "ntuples"
        }
    elif analysis == "WZDecemberAnalysis":
        selection_map = { "ntuples" : "ntuples",
                "loosepreselection" : "ntuples",
                "preselection" : "ntuples",
                "Mass3l" : "preselection",
                "Zselection" : "preselection",
                "Wselection" : "Zselection"
        }
    selection = selection.replace(";",",")
    first_selection = selection.split(",")[0].strip()
    if first_selection not in selection_map.keys():
        if "preselection" in first_selection:
            first_selection = "preselection"
        else:
            raise ValueError("Invalid selection '%s'. Valid selections are:"
                "%s" % (first_selection, selection_map.keys()))
    return selection_map[first_selection]
def getInputFilesPath(sample_name, selection, analysis, manager_path=""):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/AnalysisDatasetManager/FileInfo" % manager_path
    input_file_name = "/".join([data_path, analysis, "%s.json" %
        selection])
    input_files = UserInput.readJson(input_file_name)
    if sample_name not in input_files.keys():
        raise ValueError("Invalid input file %s. Input file must correspond"
               " to a definition in %s" % (sample_name, input_file_name))
    filename = input_files[sample_name]['file_path']
    return filename
def getCutsJsonName(selection, analysis):
    return "/".join(["Cuts", analysis, selection + ".json"]) 
def getTriggerName(sample_name, analysis, selection):
    trigger_names = ["MuonEG", "DoubleMuon", "DoubleEG", "SingleMuon", "SingleElectron"]
    if "data" in sample_name and getPreviousStep(selection, analysis) == "ntuples":
        for name in trigger_names:
            if name in sample_name:
                return "-t " + name
    return "-t MonteCarlo"
