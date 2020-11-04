import os
import argparse
import subprocess
import pandas as pd
from itertools import repeat
from multiprocessing import Pool, freeze_support

#Run Spades on a directory
def RunBileDirectory(inputDir, outputDir):
    genomeFileList = []
    gffOutList = []
    faaOutList = []
    hmmOutList = []
    for subdir, dirs, files in os.walk(inputDir):
        genomeFile = ""
        gffOut = ""
        faaOut = ""
        hmmOut = ""
        for file in files:
            if file.endswith(".fna"):
                genomeFile = os.path.join(subdir, file)
                genomeFileList.append(genomeFile)
                #outFile
                gffOut = os.path.join(outputDir, file, ".gff")
                faaOut = os.path.join(outputDir, file, ".faa")
                hmmOut = os.path.join(outputDir, file, ".tbl")
                #outFileList
                gffOutList.append(gffOut)
                faaOutList.append(faaOut)
                hmmOutList.append(hmmOut)
    RunProdigalParallel(genomeFileList, gffOutList, faaOutList, njobs)
    RunHMMerParallel(faaOutList, hmmFile, hmmOutList, ncpus, njobs)

def RunProdigalParallel(genomeFileList, gffOutList, faaOutList, njobs):
    pool = Pool(processes=njobs)
    pool.starmap(RunProdigal, zip(genomeFileList, gffOutList,  faaOutList))
    pool.close()
    pool.join()
    pool.terminate()

def RunProdigal(genomeFile, gffOut, faaOut):
    cmd = "prodigal -i " + genomeFile + " -o " + gffOut + " -a " + faaOut + " -f gff"
    subprocess.call(cmd, shell=True)

def RunHMMerParallel(fastaFileList, hmmFile, hmmOutList, ncpus, njobs):
    pool = Pool(processes=njobs)
    pool.starmap(RunHMMer, zip(fastaFileList, repeat(hmmFile), hmmOutList, repeat(ncpus)))
    pool.close()
    pool.join()
    pool.terminate()

def RunHMMer(fastaFile, hmmFile, hmmOut, ncpus):
    cmd = "hmmsearch --tblout " + hmmOut + " " + fastaFile + " " + hmmFile + " --cpu " + str(ncpus)
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Phage Assembly & Annotation')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
args = parser.parse_args()

inputDir = os.path.abspath(args.fileDir)
outputDir = os.path.abspath(args.OpDir)
njobs = int(args.jobs)
ncpus = int(args.threads)

hmmFile = "/mnt/d/Lab/Bile-Acid/phmm/CBAH.hmm"
RunBileDirectory(inputDir, outputDir)