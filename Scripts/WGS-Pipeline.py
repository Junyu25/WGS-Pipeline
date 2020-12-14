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

#Run on outDir's Spades assemble out put    
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




#Run Spades in parallel
def RunSpadesParallel(R1List, R2List, outFileList, jobs, threads):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunSpades, zip(R1List, R2List, outFileList, repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

#SPAdes Assembling
def RunSpades(R1, R2, OutDir, threads):
    os.makedirs(OutDir, 0o777, True)
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "spades.py --meta -1 " + R1 + " -2 " + R2 + " -o " + OutDir + " -t " + str(threads)
    subprocess.call(cmd, shell=True)


#Run Bandage in parallel
def RunBandageParallel(fileList, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=4)
    pool.starmap(RunBandage, zip(fileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()
#Bandage image CD1382_FDSW202399938-1r/assembly_graph.fastg CD1382.jpg
#Bandage Preview
def RunBandage(InFile, OutFile):
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "Bandage image " + os.path.join(InFile, "assemble", "assembly_graph.fastg") + " " + OutFile
    subprocess.call(cmd, shell=True)

#Run Prokka in parallel
def RunProkkaParallel(fileList, outFileList, prefixList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=4)
    pool.starmap(RunProkka, zip(fileList, outFileList, prefixList))
    pool.close()
    pool.join()
    pool.terminate()

def RunProkka(fasta, outDir, prefix):
    cmd = "prokka --kingdom Viruses --genus viral --hmms /home/junyuchen/Lab/Phage-SOP/Database/VOGDB/VOGDB_m.hmm --addgenes --prefix " + prefix + " --outdir " + outDir + " --force " + fasta
    print(cmd)
    subprocess.call(cmd, shell=True)




parser = argparse.ArgumentParser(description='WGS Assembly & Annotation')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-d', '--database', dest='Database', type=str, required=False, default='/home/malab/databases_of_malab/nr/nr',
                    help="the nr database path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
parser.add_argument('-l', '--length', dest='length', type=str, required=False, default='10000',
                    help="the length to filter contigs")
parser.add_argument('-F', '--sepF', dest='sp1', type=str, required=False, default='_1.clean.fq.gz',
                    help="It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.")
parser.add_argument('-R', '--sepR', dest='sp2', type=str, required=False, default='_2.clean.fq.gz',
                    help="It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.")
args = parser.parse_args()

inputDir = os.path.abspath(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)
nr = os.path.abspath(args.Database)
jobs = int(args.jobs)
threads = int(args.threads)
#definate length of a phage genome
dlen = int(args.length)



df = pd.read_table(inputDir)
outFileList = df["SampleID"].tolist()
R1List = df["forward-absolute-filepath"].tolist()
R2List = df["reverse-absolute-filepath"].tolist()