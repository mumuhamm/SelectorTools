#!/usr/bin/env python

import argparse
import makeFileList 
from python import UserInput,OutputTools
from python import ConfigureJobs
import datetime
import logging
import os
import shutil
import glob
import tarfile
import math
import re
import subprocess

def getComLineArgs():
    parser = UserInput.getDefaultParser(False)
    parser.add_argument("-d", "--submit_dir", type=str,
                        required=True, help="Output directory")
    parser.add_argument("-n", "--files_per_job", type=int,
                        default=3, help="Number of files per job")
    parser.add_argument("--nosubmit", action='store_true',
                        help="Just make the directory, don't submit")
    parser.add_argument("--local", action='store_true',
                        help="Use local files, e.g., file_path")
    parser.add_argument("--input_tier", type=str,
        default="", help="Selection stage of input files")
    parser.add_argument("-q", "--queue", type=str,
        default="longlunch", help="lxplus queue, or 'uw' for wisconsin settings")
    parser.add_argument("--force", action='store_true',
        help="Force overwrite of existing directories")
    return vars(parser.parse_args())

def makeSubmitDir(submit_dir, force):
    log_dir = submit_dir + "/logs"
    if os.path.isdir(submit_dir) and force:
        logging.warning("Overwriting directory %s" % submit_dir)
        shutil.rmtree(submit_dir)
    elif os.path.isdir(submit_dir):
       raise IOError("Submit directory %s already exists! Use --force to overrite." % submit_dir) 
    os.makedirs(log_dir)

def copyLibs():
    libdir = "lib"
    if os.path.isdir(libdir):
        shutil.rmtree(libdir)
    os.mkdir(libdir)
    
    cmssw_libdir = "/".join([os.environ["CMSSW_BASE"], libdir, os.environ["SCRAM_ARCH"], "*VVAnalysis*"])
    for i in glob.glob(cmssw_libdir):
        shutil.copyfile(i, '/'.join([libdir, os.path.basename(i)]))

# Needed on lxplus and uwlogin, where the afs permissions are set
# very tight and don't let condor access some dumb file it needs
# TODO: Understand what it needs and why
def modifyAFSPermissions():
    if not re.findall("system:anyuser *.*l", subprocess.check_output(["fs", "la"])):
        subprocess.call(["find", os.environ["CMSSW_BASE"], "-type", "d",
            "-exec", "fs", "setacl", "-dir", "{}", "-acl", "system:anyuser", "rl", ";"])
        raise OSError("AFS permissions have been relaxed for condor submission. You should recompile and resubmit")

def copyDatasetManagerFiles(analysis):
    manager_name = "AnalysisDatasetManager"
    manager_path = ConfigureJobs.getManagerPath()

    if os.path.isdir(manager_name):
        shutil.rmtree(manager_name)

    info_dir = manager_name+"/FileInfo"
    plot_dir = manager_name+"/PlotObjects"
    os.makedirs(info_dir)
    os.makedirs(plot_dir)

    for paths in [[manager_path, manager_name, "Utilities"],
                    [manager_path, info_dir, "data"],
                    [manager_path, info_dir, "montecarlo"],
                    [manager_path, plot_dir, analysis]]:
        path = '/'.join(paths)
        if os.path.isdir(path):
            shutil.copytree(path, '/'.join(paths[1:]))
        else:
            for d in glob.glob(path+'*'):
                shutil.copy(d, d.replace(manager_path+"/", ""))

# TODO: Check if this is needed at UW. I think it isn't
def copyGridCertificate():
    proxypath = "/tmp/x509up_u%s" % os.getuid() if not \
                ("X509_USER_PROXY" in os.environ and os.path.isfile(os.environ["X509_USER_PROXY"])) \
            else os.environ["X509_USER_PROXY"] 
    shutil.copy(proxypath, "userproxy")

def tarAnalysisInfo(condor_dir, tarball_name):
    tarname = condor_dir+"/"+tarball_name
    with tarfile.open(tarname, "w:gz") as tar:
        tar.add("Utilities")
        tar.add("data")
        tar.add("lib")
        tar.add("AnalysisDatasetManager")
    shutil.rmtree("lib")
    shutil.rmtree("AnalysisDatasetManager")

def getUWCondorSettings():
    return """# (wisconsin-specific) tell glideins to run job with access to cvmfs (via parrot)
        +RequiresCVMFS       = True
        +RequiresSharedFS    = True
        +IsFastQueueJob      = True
        Requirements         = TARGET.Arch == "X86_64" && IsSlowSlot=!=true && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.HasParrotCVMFS=?=true || (TARGET.UWCMS_CVMFS_Exists  && TARGET.CMS_CVMFS_Exists))
    """

def writeSubmitFile(submit_dir, analysis, selection, input_tier, queue, filelist, numfiles, nPerJob, selArgs):
    template_dict = {
        "analysis" : analysis,
        "selection" : selection,
        "input_tier" : input_tier,
        "queue" : ('+JobFlavour = "%s"\n+AccountingGroup = "group_u_CMST3.all"' % queue) \
                if queue != 'uw' else getUWCondorSettings(),
        "filelist" : filelist.split(".txt")[0],
        "nPerJob" : nPerJob,
        "nJobs" : int(math.ceil(numfiles/nPerJob)),
        "extraArgs" : "--debug" if not selArgs else ("--debug --selectorArgs %s" % " ".join(selArgs))
    }

    template = "Templates/CondorSubmit/submit_template.jdl"
    outfile = "/".join([submit_dir, "submit.jdl"])
    ConfigureJobs.fillTemplatedFile(template, outfile, template_dict)

def writeWrapperFile(submit_dir, tarball_name):
    template_dict = { 
            "CMSSW_RELEASE_BASE" : os.environ["CMSSW_RELEASE_BASE"],
            "tarball" : tarball_name
    }
    template = "Templates/CondorSubmit/wrapRunSelector.sh"
    outfile = "/".join([submit_dir, "wrapRunSelector.sh"])
    ConfigureJobs.fillTemplatedFile(template, outfile, template_dict)

def writeMetaInfo(submit_dir, filename):
    with open("/".join([submit_dir, filename]), "w") as metafile:
        metafile.write(OutputTools.getScriptCall()+"\n\n")
        metafile.write("Script called at %s \n" % datetime.datetime.now())
        metafile.write("git hash: " + OutputTools.gitHash()+"\n")
        metafile.write("git diff: " + OutputTools.gitDiff()+"\n")

def submitDASFilesToCondor(filenames, submit_dir, analysis, selection, input_tier, queue, numPerJob, force, das, selArgs):
    makeSubmitDir(submit_dir, force)
    copyLibs()
    copyDatasetManagerFiles(analysis)
    modifyAFSPermissions()

    filelist_name = '_'.join(filenames[:max(len(filenames), 4)])
    filelist_name = filelist_name.replace("*", "ALL")
    filelist = '/'.join([submit_dir, filelist_name+'_filelist.txt'])
    numfiles = makeFileList.makeFileList(filenames, filelist, analysis, input_tier, das)
    writeSubmitFile(submit_dir, analysis, selection, input_tier, queue, filelist_name, numfiles, numPerJob, selArgs)

    tarball_name = '_'.join([analysis, "AnalysisCode.tgz"])
    writeWrapperFile(submit_dir, tarball_name)
    tarAnalysisInfo(submit_dir, tarball_name)
    writeMetaInfo(submit_dir, "metaInfo.txt")

def main():
    args = getComLineArgs()
    submitDASFilesToCondor(args['filenames'], args['submit_dir'], args['analysis'], 
        args['selection'], args['input_tier'], args['queue'], args['files_per_job'], args['force'], 
        not args['local'], args['selectorArgs'])

if __name__ == "__main__":
    main()
