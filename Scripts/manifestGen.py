import os
import argparse
import pandas as pd

#Run Spades on a directory
def RunSpadesDirectory(inputDir):
    path = pd.DataFrame()
    for subdir, dirs, files in os.walk(inputDir):
        R1 = ""
        R2 = ""
        for file in files:
            if file.endswith(r1_end):
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file[:-len(r1_end)]+r2_end)
                sampleStr = file.replace(r1_end, "")
    path = path.append({'SampleID':str(sampleStr), "forward-absolute-filepath":str(R1), "reverse-absolute-filepath":str(R2)}, ignore_index=True)

parser = argparse.ArgumentParser(description='Manifest Gen')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=False, default='./',
                    help="the output path of reads")
parser.add_argument('-F', '--sepF', dest='sp1', type=str, required=False, default='_1.clean.fq.gz',
                    help="It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.")
parser.add_argument('-R', '--sepR', dest='sp2', type=str, required=False, default='_2.clean.fq.gz',
                    help="It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.")
args = parser.parse_args()


r1_end = args.sp1
r2_end = args.sp2
inputDir = os.path.abspath(args.fileDir)
OpDir = os.path.abspath(args.OpDir)
path = RunSpadesDirectory(inputDir)
path.to_csv(OpDir+"/PathTable.tsv", sep="\t", index=False, encoding = "utf-8")
