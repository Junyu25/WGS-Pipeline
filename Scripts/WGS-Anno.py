'''
Copyright {2020} Junyu Chen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import os
import argparse
import subprocess
import pandas as pd
#from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support

# Don't run in the Dir 
# It's better to run in the manifest profile.

## Fastp pair end tirmming and merging
def RunFastp(R1, R2, prefix, OutDir, threads):
    cmd = "fastp --in1 " + R1 + " --in2 " + R2 + " --out1 " + os.path.join(OutDir, prefix + "_R1.fastq") + " --out2 " + os.path.join(fastpDir, prefix + "_R2.fastq") + \
    " --thread " + str(threads) + \
    " --html " + os.path.join(OutDir, prefix + ".html") + " --json " + os.path.join(OutDir, prefix + ".json") + " --report_title " + prefix + "-fastq-merge-report"
    subprocess.call(cmd, shell=True)
## Run fastp in parallel
def RunFastpParallel(R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunFastp, zip(R1List, R2List, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

## RunUnicycler
def RunUnicycler(R1, R2, prefix, OutDir, threads):
    cmd = "unicycler -1 " + R1 + " -2 " + R2 + " -o " + os.path.join(OutDir, prefix) + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
def RunUnicyclerParallel(R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunUnicycler, zip(R1List, R2List, prefixList, repeat(OutDir), repeat(threads)))
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
    cmd = "quast " + fasta + " -1 " + R1 + " -2 " + R2 + " -o " + os.path.join(OutDir, prefix) + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
## RunQuastParallel
def RunQuastParallel(fastaList, R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunQuast, zip(fastaList, R1List, R2List, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

def RunBandage(InFile, prefix, OutDir):
    #InFile: A graph file of any type supported by Bandage
    #OutFile: The image file to be created (must end in '.jpg', '.png' or '.svg')
    cmd = "Bandage image " + InFile + " " + os.path.join(OutDir, prefix + "_graph.png")
    subprocess.call(cmd, shell=True)
#Run Bandage in parallel
def RunBandageParallel(fileList, prefixList, OutDir, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunBandage, zip(fileList, prefixList, repeat(OutDir)))
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

    
def parseRaw(prefixList, fastpDir):
    R1_cleanList = []
    R2_cleanList = []
    for prefix in prefixList:
        R1_cleanList.append(os.path.join(fastpDir, prefix + "_R1.fastq"))
        R2_cleanList.append(os.path.join(fastpDir, prefix + "_R2.fastq"))
    return R1_cleanList, R2_cleanList

def parseAssembly(prefixList, UnicyclerDir):
    assembleList = []
    graphList = []
    for prefix in prefixList:
        assembleList.append(os.path.join(UnicyclerDir, prefix, "assembly.fasta"))
        graphList.append(os.path.join(UnicyclerDir, prefix, "assembly.gfa")) #assembly.gfa
    #print(assembleList)
    #print(graphList)
    return assembleList, graphList
'''
.fna	Nucleotide FASTA file of the input contig sequences.
.faa	Protein FASTA file of the translated CDS sequences.
.ffn	Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
'''
def parseProkka(prefixList, ProkkaDir):
    fnaList = []
    faaList = []
    ffnList = []
    for prefix in prefixList:
        fnaList.append(os.path.join(ProkkaDir, prefix, prefix + ".fna"))
        faaList.append(os.path.join(ProkkaDir, prefix, prefix + ".faa"))
        ffnList.append(os.path.join(ProkkaDir, prefix, prefix + ".ffn"))
    return fnaList, faaList, ffnList

## Run kraken2
def RunKraken2(fasta, database, prefix, OutDir, threads):
    cmd = "kraken2 -db " + database + " " + fasta + " --report " + os.path.join(OutDir, prefix + "_kraken2.tsv") + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
#Run kraken2 in parallel
def RunKraken2Parallel(fileList, database, prefixList, OutDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunKraken2, zip(fileList, repeat(database), prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()
#parse Kraken2 result
def parseKraken2(report):
    df = pd.read_table(report)
    df.columns = ["%", "reads", "lreads", "lvl", "tid", "name"]
    df1 = df.loc[df["lvl"] == "S"]
    df2 = df1.loc[df1["%"] == df1["%"].max()]
    df2 = df2.reset_index()
    tid = df2["tid"][0]
    handle = Entrez.efetch(db="Taxonomy", id=str(tid) , retmode="xml")
    records = Entrez.read(handle)
    taxa = records[0]["Lineage"] + "; " + records[0]["ScientificName"]
    return tid, taxa

## Run ABRicate
def RunABRicate(fasta, prefix, OutDir, threads):
    cmd = "abricate " + fasta + " > " + os.path.join(OutDir, prefix + "_abricate.txt") + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
#Run ABRicate in parallel
def RunABRicateaParallel(fileList, prefixList, OutDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunABRicate, zip(fileList, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

## Run KofamScan
def RunKofamScan(faaFile, prefix, OutDir, profile, koList, threads):
    cmd = "exec_annotation " + faaFile + " -o " + os.path.join(OutDir, prefix + "_kofam_scan.tsv") + \
        " --profile " + profile + " --ko-list " + koList + " --cpu " + str(threads) + \
        " -f detail-tsv --e-value=0.01 --tmp-dir " + os.path.join(OutDir, prefix + "_tmp") #setting tem-dir for parallel processing
    print(cmd+"\n")
    subprocess.call(cmd, shell=True)
def RunKofamScanParallel(faaList, prefixList, OutDir, profile, koList, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunKofamScan, zip(faaList, prefixList, repeat(OutDir), repeat(profile), repeat(koList), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

## Run dbCAN
def RunDbCAN(faaFile, prefix, OutDir, databseDir, threads):
    cmd = "run_dbcan.py " + faaFile + " protein --out_dir " + os.path.join(OutDir, prefix) + " --db_dir " + databseDir + \
        " --dia_cpu " + str(threads) + " --hmm_cpu " + str(threads) + " --hotpep_cpu " + str(threads) + " --tf_cpu " + str(threads)
    subprocess.call(cmd, shell=True)
def RunDbCANParallel(faaList, prefixList, OutDir, databseDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunDbCAN, zip(faaList, prefixList, repeat(OutDir), repeat(databseDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()


parser = argparse.ArgumentParser(description='WGS Assembly & Annotation')
parser.add_argument('-i', '--input', dest='InFile', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-kd', '--database', dest='kraken2_db', type=str, required=False, default='/home/junyuchen/3-Resources/Databases/k2_standard_20201202',
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
df = pd.read_csv(InFile)
prefixList = df["SampleID"].tolist()
fnaList = df["fna"].tolist()
faaList = df["faa"].tolist()
ffnList = df["ffn"].tolist()

## init out dir
if os.path.exists(OutDir) == 0:
    os.makedirs(OutDir, 0o777, True)
Kraken2Dir = os.path.join(OutDir, "kraken2")
if os.path.exists(Kraken2Dir) == 0:
    os.makedirs(Kraken2Dir, 0o777, True)
ABRicateDir = os.path.join(OutDir, "abricate")
if os.path.exists(ABRicateDir) == 0:
    os.makedirs(ABRicateDir, 0o777, True)
KofamScanDir = os.path.join(OutDir, "kofam_scan")
if os.path.exists(KofamScanDir) == 0:
    os.makedirs(KofamScanDir, 0o777, True)
dbCANDir = os.path.join(OutDir, "dbcan")
if os.path.exists(dbCANDir) == 0:
    os.makedirs(dbCANDir, 0o777, True)
BileAcidDir = os.path.join(OutDir, "bile_acid")
if os.path.exists(BileAcidDir) == 0:
    os.makedirs(BileAcidDir, 0o777, True)


RunKraken2Parallel(fnaList, kraken2_db, prefixList, Kraken2Dir, threads, jobs)

RunABRicateaParallel(fnaList, prefixList, ABRicateDir, threads, jobs)

profile1 = "/home/junyuchen/1-Projects/WGS-Pipeline/database/bile_acid/profiles"
koList1 = "/home/junyuchen/1-Projects/WGS-Pipeline/database/bile_acid/bile_acid_ko_list"
RunKofamScanParallel(faaList, prefixList, BileAcidDir, profile1, koList1, threads, jobs)

profile = "/home/junyuchen/1-Projects/WGS-Pipeline/database/kofam_scan/profiles"
koList = "/home/junyuchen/1-Projects/WGS-Pipeline/database/kofam_scan/ko_list"
RunKofamScanParallel(faaList, prefixList, KofamScanDir, profile, koList, threads, jobs)

databseDir = "/home/junyuchen/1-Projects/WGS-Pipeline/database/dbCAN2_db"
RunDbCANParallel(faaList, prefixList, dbCANDir, databseDir, threads, jobs)