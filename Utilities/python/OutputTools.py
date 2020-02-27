import ROOT
import datetime
import subprocess
import sys
import logging

def addMetaInfo(fOut):
    metaInfo = fOut.mkdir("MetaInfo")
    metaInfo.cd()
    time = ROOT.TNamed("datetime", str(datetime.datetime.now()))
    command = ROOT.TNamed("command", getScriptCall())
    try: 
        githash = ROOT.TNamed("githash", gitHash())
        gitdiff = ROOT.TNamed("gitdiff", gitDiff())
        githash.Write()
        gitdiff.Write()
    except subprocess.CalledProcessError:
        logging.warning("Not a git repo, skipping git hash/diff")
        pass

    time.Write()
    command.Write()

def gitHash():
    return subprocess.check_output(['git', 'log', '-1', '--format="%H"'])

def gitDiff():
    return subprocess.check_output(['git', 'diff',])

def getScriptCall():
    return ' '.join(sys.argv)

def writeOutputListItem(item, directory):
    if item.ClassName() == "TList":
        d = directory.Get(item.GetName())
        if not d:
            d = directory.mkdir(item.GetName())
            ROOT.SetOwnership(d, False)
        for subItem in item:
            writeOutputListItem(subItem, d)
    elif hasattr(item, 'Write'):
        directory.cd()
        item.Write()
    else:
        print "Couldn't write output item:"
        print repr(item)
    directory.cd()
