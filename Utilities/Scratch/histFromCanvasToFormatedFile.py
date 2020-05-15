import argparse
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_files", nargs="*", type=str, required=True)
parser.add_argument("-o", "--output_file", type=str, required=True)
parser.add_argument("-p", "--process_names", nargs="*", type=str, required=True)
parser.add_argument("-c", "--canvas", type=str, default="canvas")
parser.add_argument("-b", "--hist_name", type=str, required=True)
args = parser.parse_args()

outfile = ROOT.TFile(args.output_file, "RECREATE")

for f, process in zip(args.input_files, args.process_names):
    rtfile = ROOT.TFile.Open(f)
    hist = 0
    outfile.mkdir(process)
    outfile.cd(process)
    canvas = rtfile.Get(args.canvas)
    if not canvas:
        raise ValueError("Could not find canvas %s in file" % args.canvas)
    for h in canvas.GetListOfPrimitives():
        if h.InheritsFrom("TH1"):
            hist = h.Clone(args.hist_name)
            break
    hist.Write()
    outfile.cd()


