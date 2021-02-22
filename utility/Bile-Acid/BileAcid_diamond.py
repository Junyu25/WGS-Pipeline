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
    blastpOutList = []
    for subdir, dirs, files in os.walk(inputDir):
        genomeFile = ""
        gffOut = ""
        faaOut = ""
        blastpOut = ""
        for file in files:
            if file.endswith(".fna"):
                genomeFile = os.path.join(subdir, file)
                genomeFileList.append(genomeFile)
                #outFile
                gffOut = os.path.join(outputDir, file + ".gff")
                faaOut = os.path.join(outputDir, file + ".faa")
                blastpOut = os.path.join(outputDir, file + ".tsv")
                #outFileList
                gffOutList.append(gffOut)
                faaOutList.append(faaOut)
                blastpOutList.append(blastpOut)
    RunProdigalParallel(genomeFileList, gffOutList, faaOutList, njobs)
    RunDiamondParallel(faaOutList, db, njobs, nthreads, blastpOutList)

def RunProdigalParallel(genomeFileList, gffOutList, faaOutList, njobs):
    pool = Pool(processes=njobs)
    pool.starmap(RunProdigal, zip(genomeFileList, gffOutList,  faaOutList))
    pool.close()
    pool.join()
    pool.terminate()

def RunProdigal(genomeFile, gffOut, faaOut):
    cmd = "prodigal -i " + genomeFile + " -o " + gffOut + " -a " + faaOut + " -f gff"
    subprocess.call(cmd, shell=True)

def RunDiamondParallel(fastaList, db, jobs, threads, outFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunDiamond, zip(fastaList, repeat(db), repeat(threads), outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunDiamond(fasta, db, threads, OutFile):
    cmd = "diamond blastp -q " + fasta  + " -o " + OutFile + " --evalue 1.0 --max-target-seqs 1 --outfmt 6 --db " + db + " -p " + str(threads) 
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
parser.add_argument('-d', '--databse', dest='databse', type=str, required=False, default='/mnt/d/Lab/Bile-Acid/db/HSDH',
                    help="the number of threads run for a job")
args = parser.parse_args()

inputDir = os.path.abspath(args.fileDir)
outputDir = os.path.abspath(args.OpDir)
njobs = int(args.jobs)
nthreads = int(args.threads)
db = os.path.abspath(args.databse)

RunBileDirectory(inputDir, outputDir)