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
#import argparse
import subprocess
#import pandas as pd
#from Bio import SeqIO
#from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support


def RunDiamond(fasta, prefix, dbDir, OutDir, threads):
    if os.path.exists(OutDir) == 0:
        os.makedirs(OutDir, 0o777, True)
    cmd = "diamond blastp --db " + dbDir + " --query " + fasta + " --out " + os.path.join(OutDir, prefix + "_blasp.tsv") + " --evalue 1e-05 --outfmt 6 --max-target-seqs 1" + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
#Run Diamond in parallel
def RunDiamondParallel(fileList, prefixList, dbDir, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunDiamond, zip(fileList, prefixList, repeat(dbDir), repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

## Run Prokka
def RunProkka(fasta, prefix, OutDir, threads):
    cmd = "prokka --addgenes --prefix " + prefix + " --outdir " + os.path.join(OutDir, prefix) + \
        " --force " + fasta + " --cpus " + str(threads)
    #print(cmd)
    subprocess.call(cmd, shell=True)
## Run Prokka in parallel
def RunProkkaParallel(fileList, prefixList, OutDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunProkka, zip(fileList, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

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

## Run ABRicate
def RunABRicate(fasta, prefix, OutDir, threads):
    cmd = "abricate " + fasta + " > " + os.path.join(OutDir, prefix + "_abricate.txt") + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)
## Run ABRicate in parallel
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
#def parseDbCAN():
