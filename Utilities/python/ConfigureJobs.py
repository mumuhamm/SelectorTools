import datetime
import UserInput
import fnmatch 
import glob
import subprocess
import os
import json
import array
import string
import socket
import logging
import logging
#try:
import configparser
#except:
    #import ConfigParser as configparser
    #from six.moves import configparser

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

def getBinning(variable='MTWZ', isVBS=True, isHiggs=False):
    if variable == 'MTWZ':
        if isVBS:
            if isHiggs:
                return [0,50,100,150,200,250,300,400,500,700,1000,1500,2000]
            return [0,100,200,300,400,500,700,1000,1500,2000]
        return [0,50,100,200,300,400,500,700,1000,1200]
    return []

def getChannels(analysis='WZ'):
    if analysis == 'WZ':
        return ["eee", "eem", "emm", "mmm"]

def getManagerName():
    config_name = ""
    try:
        config_name = "Templates/config.%s" % os.getlogin()
    except OSError:
        pass
    default_name = "AnalysisDatasetManager"
    if not os.path.isfile(config_name):
        logging.warning("dataset_manager_name not specified in config file %s" % config_name)
        logging.warning("Using default '%s'" % default_name)
        return default_name
    config = configparser.ConfigParser()
    config.read_file(open(config_name))
    if "dataset_manager_name" not in config['Setup']:
        logging.warning("dataset_manager_name not specified in config file %s" % config_name)
        logging.warning("Using default '%s'" % default_name)
        return default_name
    return config['Setup']['dataset_manager_name']

def getManagerPath():
    config_name = ""
    try:
        config_name = "Templates/config.%s" % os.getlogin()
    except OSError:
        pass
    if not os.path.isfile(config_name):
        if os.path.isdir(getManagerName()):
            return '.'
        else:
            raise IOError("Failed to find valid config file. Looking for %s" 
                    % config_name)
    config = configparser.ConfigParser()
    config.read_file(open(config_name))
    if "dataset_manager_path" not in config['Setup']:
        raise ValueError("dataset_manager_path not specified in config file %s"
                        % config_name)
    return os.path.expanduser(config['Setup']['dataset_manager_path']) + "/"

def getCombinePath():
    config = configparser.ConfigParser()
    config.read_file(open("Templates/config.%s" % os.environ["USER"]))
    if "combine_path" not in config['Setup']:
        raise ValueError("dataset_manager_path not specified in config file Template/config.%s" 
                            % os.environ["USER"])
    return config['Setup']['combine_path'] + "/"
def getListOfEWKFilenames(analysis=""):
    if "ZZ4l" in analysis:
        return [
            "zz4l-powheg",
            "ggZZ4e",
            "ggZZ4m",
            "ggZZ4t",
            "ggZZ2e2mu",
            "ggZZ2e2tau",
            "ggZZ2mu2tau",
        ]
    # TODO: This is obviously WZ specific and should be updated
    return [
    #    "wz3lnu-powheg",
    # Use jet binned WZ samples for subtraction by default
        "wz3lnu-mgmlm-0j",
        "wz3lnu-mgmlm-1j",
        "wz3lnu-mgmlm-2j",
        "wz3lnu-mgmlm-3j",
        "wlljj-ewk",
        "zz4l-powheg",
        "zz4ljj-ewk",
        "zz2l2vjj-ewk",
        "tzq",
        "ttz",
        "ttw",
        "zzz",
        "wwz",
        "www",
        "ww",
        "zg",
        "ggZZ4e",
        "ggZZ4m",
        "ggZZ2e2mu",
    ]
def getListOfNonpromptFilenames():
    return ["tt-lep",
        "st-schan",
        "st-tchan-t",
        "st-tchan-tbar",
        "st-tw",
        "st-tbarw",
        #"DYm50",
        "DYm50-1j",
        "DYm50-2j",
        "DYm50-3j",
        "DYm50-4j",
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
    p = subprocess.Popen(["hdfs", "dfs", "-ls", "-h", file_path.replace("/hdfs", "")],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    out,err = p.communicate()
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
        split = line.rsplit(" ", 1)
        if len(split) != 2:
            continue
        elif "root" in split[1]:
            files.append(split[1])
    return files

# TODO: Would be good to switch the order of the last two arguments
# completely deprecate manager_path without breaking things
def getListOfFiles(filelist, selection, manager_path="", analysis=""):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/%s/FileInfo" % (manager_path, getManagerName())
    group_path = "%s/AnalysisDatasetManager/PlotGroups" % manager_path
    data_info = UserInput.readAllInfo("/".join([data_path, "data/*"]))
    mc_info = UserInput.readAllInfo("/".join([data_path, "montecarlo/*"]))
    analysis_info = UserInput.readInfo("/".join([data_path, analysis, selection])) \
        if analysis != "" else []
    valid_names = (data_info.keys() + mc_info.keys()) if not analysis_info else analysis_info.keys()
    group_names = UserInput.readAllInfo("%s/%s.py" %(group_path, analysis)) if analysis else dict()
    names = []
    for name in filelist:
        if ".root" in name:
            names.append(name)
        # Allow negative contributions
        elif name[0] == "-" :
            names.append(name)
        elif "WZxsec2016" in name:
            dataset_file = manager_path + \
                "%s/FileInfo/WZxsec2016/%s.json" % (getManagerPath(), selection)
            allnames = json.load(open(dataset_file)).keys()
            if "nodata" in name:
                nodata = [x for x in allnames if "data" not in x]
                names += nodata
            elif "data" in name:
                names += [x for x in allnames if "data" in x]
            else:
                names += allnames
        elif name in group_names:
            names += group_names[name]['Members']
        elif "*" in name:
            anti = "NOT" in name[:3]
            name = name.replace("NOT", "")
            matching = fnmatch.filter(valid_names, name)
            if anti:
                names += filter(lambda x: x not in matching, valid_names)       
            else:
                names += matching
        elif name not in valid_names and name.split("__")[0] not in valid_names:
            logging.warning("%s is not a valid name" % name)
            logging.warning("Valid names must be defined in AnalysisDatasetManager/FileInfo/(data/montecarlo)*")
            logging.debug("Vaid names are %s" % valid_names)
            continue
        else:
            names += [name]
    if not names or len(filter(lambda x: x != '', names)) == 0:
        raise RuntimeError("No processes found matching pattern '%s'" % filelist)
    return [str(i) for i in names]

def getXrdRedirector(filepath=None):
    if "eos/cms" in filepath:
        return "eoscms.cern.ch"
    elif "eos/user" in filepath:
        return "eosuser.cern.ch"

    usbased = ["wisc.edu"]
    usredir = 'cmsxrootd.fnal.gov'
    globalredir = 'cms-xrd-global.cern.ch'
    # Cluster machines may not have this env variable
    if any(i in socket.gethostname() for i in usbased):
        return usredir
    return globalredir

def fillTemplatedFile(template_file_name, out_file_name, template_dict):
    with open(template_file_name, "r") as templateFile:
        source = string.Template(templateFile.read())
        result = source.substitute(template_dict)
    with open(out_file_name, "w") as outFile:
        outFile.write(result)

def getListOfFilesWithXSec(filelist, manager_path="", selection="ntuples"):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/%s/FileInfo" % (manager_path, getManagerName())
    files = getListOfFiles(filelist, selection, manager_path)
    mc_info = UserInput.readAllInfo("/".join([data_path, "montecarlo/*"]))
    info = {}
    for file_name in files:
        if "data" in file_name.lower() or "nonprompt" in file_name.lower():
            info.update({file_name : 1})
        else:
            label = file_name.split("__")[0]
            if label[0] == "-":
                label = label[1:]
            file_info = mc_info[label]
            kfac = file_info["kfactor"] if "kfactor" in file_info.keys() else 1
            # Intended to use for, e.g., nonprompt
            if file_name[0] == "-":
                kfac *= -1
            info.update({file_name : file_info["cross_section"]*kfac})
    return info

def getListOfFilesWithPath(filelist, analysis, selection, das=True, manager_path=""):
    if manager_path is "":
        manager_path = getManagerPath()
    data_path = "%s/%s/FileInfo" % (manager_path, getManagerName())
    files = getListOfFiles(filelist, selection, manager_path, analysis)
    selection_info = UserInput.readInfo("/".join([data_path, analysis, selection]))
    info = {}
    for file_name in files:
        if das and "DAS" not in selection_info[file_name].keys():
            logging.error("DAS path not defined for file %s in analysis %s/%s" % (file_name, analysis, selection))
            continue
        elif not das and "file_path" not in selection_info[file_name].keys():
            logging.error("File_path not defined for file %s in analysis %s/%s" % (file_name, analysis, selection))
            continue
        info.update({file_name : selection_info[file_name]["DAS" if das else "file_path"]})
    return info

def getPreviousStep(selection, analysis):
    selection_map = {}
    if analysis == "WZxsec2016":
        selection_map = { "ntuples" : "ntuples",
                "loosepreselection" : "ntuples",
                "preselection" : "ntuples",
                "3LooseLeptons" : "ntuples",
                "3LooseLeptonsNoIP" : "ntuples",
                "3LooseLeptonsNoVeto" : "ntuples",
                "3TightLeptonsNoVeto" : "ntuples",
                "3MediumLeptonsNoVeto" : "ntuples",
                "preselectionLooseVeto" : "ntuples",
                "preselectionNoVeto" : "ntuples",
                "LepVetoAnd3lmass" : "preselection",
                "3lmass" : "preselection",
                "Zselection" : "3lmass",
                "Wselection" : "Zselection",
                "3lDYControl" : "Zselection",
                "3lTTbarControl" : "3lmass",
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

def getConfigFileName(config_file_name):
    for extension in ["json", "py"]:
        if os.path.isfile(".".join([config_file_name, extension])):
            return ".".join([config_file_name, extension])
    raise ValueError("Invalid configuration file. Tried to read %s which does not exist" % \
            config_file_name)

def getInputFilesPath(sample_name, selection, analysis, manager_path=""):
    if manager_path is "":
        manager_path = getManagerPath()
    if ".root" in sample_name:
        logging.info("Using simple file %s" % sample_name)
        return sample_name
    data_path = "%s/%s/FileInfo" % (manager_path, getManagerName())
    input_file_base_name = "/".join([data_path, analysis, selection])
    input_file_name = getConfigFileName(input_file_base_name)
    input_files = UserInput.readInfo(input_file_name)
    if sample_name not in input_files.keys():
        raise ValueError("Invalid input file %s. Input file must correspond"
               " to a definition in %s" % (sample_name, input_file_name))
    filename = input_files[sample_name]['file_path']
    return filename

def getCutsJsonName(selection, analysis):
    return "/".join(["Cuts", analysis, selection])

def getTriggerName(sample_name, analysis, selection):
    trigger_names = ["MuonEG", "DoubleMuon", "DoubleEG", "SingleMuon", "SingleElectron"]
    if "data" in sample_name and getPreviousStep(selection, analysis) == "ntuples":
        for name in trigger_names:
            if name in sample_name:
                return "-t " + name
    return "-t MonteCarlo"
