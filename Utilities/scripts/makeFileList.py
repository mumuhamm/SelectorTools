#!/usr/bin/env python

import argparse
import subprocess
from python import UserInput
from python import ConfigureJobs
import logging
import os
import random
import glob
import logging

def getComLineArgs():
    parser = UserInput.getDefaultParser(False)
    parser.add_argument("-o", "--output_file", type=str,
                        required=True, help="Name of output text file")
    parser.add_argument("--das", action='store_true',
                        help="Read files from DAS (default local, e.g., from file_path")
    return vars(parser.parse_args())

def getFilesWithName(name, path, das=True):
    if das:
        files = subprocess.check_output(["dasgoclient", "--query=file dataset=%s" % path])
        if files:
            files = files.split("\n")
            files = filter(lambda x: "/store" in x[:7], files)
    else:
        files = glob.glob(path)

    if files:
        files = ["@".join([name, f.strip()+'\n']) for f in files]
    return files

def makeFileList(filenames, output_file, analysis, selection, das):
    name_path_map = ConfigureJobs.getListOfFilesWithPath(filenames, analysis, selection, das)

    files = []
    for name, path in name_path_map.iteritems():
        logging.debug("Trying to add name %s with path %s" % (name, path))
        try:
            files.extend(getFilesWithName(name, path, das))
        except subprocess.CalledProcessError:
            logging.warning("Failed to find files for dataset %s with DAS path %s. Skipping!" % (name, path))
    
    # Randomize so one job doesn't tend to take longer than others
    random.shuffle(files)
    with open(output_file, "w") as outfile:
        outfile.writelines(files)
    if not len(files):
        raise RuntimeError("Failed to find any files for the given data sets")

    return len(files)

def main():
    args = getComLineArgs()
    makeFileList(args['filenames'], args['output_file'], args['analysis'], args['selection'], args['das'])

if __name__ == "__main__":
    main()

