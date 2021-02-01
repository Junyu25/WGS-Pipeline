import os
import argparse
import subprocess
import pandas as pd
#from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support


#Run Spades on a directory
def RunSpadesDirectory(inputDir, ouputDir):
    R1List = []
    R2List = []
    outFileList = []
    SpadesFileList = []
    SpadesOutList = []
    BandageOutList = []
    for subdir, dirs, files in os.walk(inputDir):
        R1 = ""
        R2 = ""
        outputFilePath = ""
        SpadesFilePath = ""
        SpadesOutDir = ""
        for file in files:
            if file.endswith(r1_end):
                outputFilePath = os.path.join(ouputDir, sampleStr)
                outFileList.append(outputFilePath)
                #SPAdes file out
                SpadesOutDir = os.path.join(ouputDir, sampleStr, "assemble")
                SpadesFilePath = os.path.join(SpadesOutDir, "scaffolds.fasta")
                #SPAdes list
                SpadesOutList.append(SpadesOutDir)
                SpadesFileList.append(SpadesFilePath)
                #Bandage
                BandageOutList.append(os.path.join(ouputDir, sampleStr, "preview.png"))
                #make out dir for every run
                os.makedirs(outputFilePath, 0o777, True)

    RunSpadesParallel(R1List, R2List, SpadesOutList, jobs, threads)
    RunBandageParallel(outFileList, BandageOutList)
    #RunProkkaParallel(SpadesFilePath, outFileList, SpadesFilePath) #prefix?

#Run on outDir's Spades assemble output    
def RunProkkaDirectory(Dir):
    prefixList = []
    contigsFileList = []
    contigsOutDirList = []
    for run in os.listdir(Dir):
        for assemble in os.listdir(os.path.join(Dir, run)):
            contigsOutDir = ""
            contigsOutPath = ""
            SpadesFilePath = os.path.join(Dir, run, assemble, "scaffolds.fasta")
            print(SpadesFilePath)
            
        RunProkkaParallel(contigsFileList, contigsOutDirList, prefixList)

# Don't run in the Dir 
# It's better to run in the manifest profile.

## Fastp pair end tirmming and merging
def RunFastp(R1, R2, prefix, OutDir, threads):
    fastpDir = os.path.join(OutDir, "fastp")
    if os.path.exists(OutDir) == 0:
        os.makedirs(OutDir, 0o777, True)
    if os.path.exists(fastpDir) == 0:
        os.makedirs(fastpDir, 0o777, True)
    cmd = "fastp --in1 " + R1 + " --in2 " + R2 + " --out1 " + os.path.join(fastpDir, prefix + "_R1.fastq") + " --out2 " + os.path.join(fastpDir, prefix + "_R2.fastq") + \
    " --thread " + str(threads) + \
    " --html " + os.path.join(fastpDir, prefix + ".html") + " --json " + os.path.join(fastpDir, prefix + ".json") + " --report_title " + prefix + "-fastq-merge-report"
    subprocess.call(cmd, shell=True)
## Run fastp in parallel
def RunFastpParallel(R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunFastp, zip(R1List, R2List, prefixList, repeat(thread), repeat(OutDir)))
    pool.close()
    pool.join()
    pool.terminate()

## RunUnicycler
def RunUnicycler(R1, R2, prefix, OutDir, threads):
    cmd = "unicycler -1 " + R1 + " -2 " + R2 + " -o " + os.path.join(OutDir, prefix) + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
def RunUnicyclerParallel(R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunUnicycler, zip(R1List, R2List, prefixList, repeat(OutDir), repeat(thread)))
    pool.close()
    pool.join()
    pool.terminate()

## SPAdes Assembling
def RunSpades(R1, R2, OutDir, threads):
    #os.makedirs(OutDir, 0o777, True)
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "spades.py --meta -1 " + R1 + " -2 " + R2 + " -o " + OutDir + " -t " + str(threads)
    subprocess.call(cmd, shell=True)
## Run Spades in parallel
def RunSpadesParallel(R1List, R2List, outFileList, jobs, threads):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunSpades, zip(R1List, R2List, outFileList, repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()


def RunQuast(fasta, R1, R2, prefix, OutDir, threads):
    cmd = "quast.py " + fasta + " -1 " + R1 + " -2 " + R2 + " -o " + os.path.join(OutDir, prefix) + " --threads " + str(threads)
## RunQuastParallel
def RunQuastParallel(fastaList, R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunQuast, zip(fastaList, R1List, R2List, prefixList, repeat(thread), repeat(OutDir)))
    pool.close()
    pool.join()
    pool.terminate()

def RunBandage(InFile, OutFile):
    #InFile: A graph file of any type supported by Bandage
    #OutFile: The image file to be created (must end in '.jpg', '.png' or '.svg')
    cmd = "Bandage image " + os.path.join(InFile, "assemble", "assembly_graph.fastg") + " " + OutFile
    subprocess.call(cmd, shell=True)
#Run Bandage in parallel
def RunBandageParallel(fileList, outFileList, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunBandage, zip(fileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

## Run Prokka
def RunProkka(fasta, prefix, OutDir, threads):
    cmd = "prokka --addgenes --prefix " + prefix + " --outdir " + os.path.join(OutDir, prefix) + " --force " + fasta + " --cpus " + str(threads)
    print(cmd)
    subprocess.call(cmd, shell=True)
#Run Prokka in parallel
def RunProkkaParallel(fileList, prefixList, OutDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunProkka, zip(fileList, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()




parser = argparse.ArgumentParser(description='WGS Assembly & Annotation')
parser.add_argument('-i', '--input', dest='InFile', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-kd', '--database', dest='kraken2_db', type=str, required=False, default='/home/malab/databases_of_malab/nr/nr',
                    help="the nr database path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
parser.add_argument('-F', '--sepF', dest='sp1', type=str, required=False, default='_1.clean.fq.gz',
                    help="It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.")
parser.add_argument('-R', '--sepR', dest='sp2', type=str, required=False, default='_2.clean.fq.gz',
                    help="It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.")
args = parser.parse_args()

InFile = os.path.abspath(args.InFile)
OutDir = os.path.abspath(args.OutDir)
kraken2_db = os.path.abspath(args.kraken2_db)
jobs = int(args.jobs)
threads = int(args.threads)

## process manifest
df = pd.read_table(InFile)
outFileList = df["SampleID"].tolist()
R1List = df["forward-absolute-filepath"].tolist()
R2List = df["reverse-absolute-filepath"].tolist()
## init out dir
UnicyclerDir = os.path.join(OutDir, "unicycler")
if os.path.exists(UnicyclerDir) == 0:
    os.makedirs(UnicyclerDir, 0o777, True)
QuastrDir = os.path.join(OutDir, "quast")
if os.path.exists(QuastrDir) == 0:
    os.makedirs(QuastrDir, 0o777, True)
BandageDir = os.path.join(OutDir, "bandage")
if os.path.exists(BandageDir) == 0:
    os.makedirs(BandageDir, 0o777, True)
Kraken2Dir = os.path.join(OutDir, "kraken2")
if os.path.exists(Kraken2Dir) == 0:
    os.makedirs(Kraken2Dir, 0o777, True)
ProkkaDir = os.path.join(OutDir, "prokka")
if os.path.exists(ProkkaDir) == 0:
    os.makedirs(ProkkaDir, 0o777, True)
ABRicateDir = os.path.join(OutDir, "abricate")
if os.path.exists(ABRicateDir) == 0:
    os.makedirs(ABRicateDir, 0o777, True)