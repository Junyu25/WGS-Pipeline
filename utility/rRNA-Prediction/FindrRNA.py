# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 09:52:16 2019

@author: Junyu
"""

import os
from BCBio import GFF
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

dir = "/home/junyuchen/Lab/16S-Prediction/PROKKA_07062020/"

limit_info = dict(gff_type = ["rRNA"])
        
def FindrRNA(path, tempFile):
    
    in_handle = open(path)

    for record in GFF.parse(in_handle, limit_info=limit_info):

        if len(record.features) != 0:
            for i in range(len(record.features)):

                if '23S ribosomal RNA' in record.features[i].qualifiers["product"][0]:
                    #print(record.features[i])
                    #print(record.id)
                    #print(record.features[i].location.strand)
                    #print(record.features[i].location.start)
                    #print(record.features[i].location.end)
                    #print(type(record.features[i].location))
                    #feature_seq = record.seq[record.features[i].location.start:record.features[i].location.end].reverse_complement()
                    #print(feature_seq)
                    '''
                    seq_23S = SeqRecord(record.features[i].location.extract(record.seq),\
                                    id=record.features[i].id,\
                                    description=str(record.id +'-'+ record.features[i].qualifiers["product"][0]))
                    '''
                    seq_23S = SeqRecord(record.features[i].location.extract(record.seq),\
                                    id=tempFile[:-4],\
                                    description=str(record.id +'-'+record.features[i].id+'-'+ record.features[i].qualifiers["product"][0]+'-'+str(len(record.seq))))
                    #print(seq_23S)
                    #SeqIO.write(seq_23S, path[:-4]+"_23S_rRNA.fasta", "fasta")
                    SeqIO.write(seq_23S, "/home/junyuchen/Lab/16S-Prediction/"+tempFile[:-4]+"_23S_rRNA.fasta", "fasta")
                elif '16S ribosomal RNA' in record.features[i].qualifiers["product"][0]:
                    #print(record.features[i])
                    #print(record.id)
                    seq_16S = SeqRecord(record.features[i].location.extract(record.seq),\
                                    id=tempFile[:-4],\
                                    description=str(record.id +'-'+record.features[i].id+'-'+ record.features[i].qualifiers["product"][0]+'-'+str(len(record.seq))))
                
                    #SeqIO.write(seq_16S, path[:-4]+"_16S_rRNA.fasta", "fasta")
                    SeqIO.write(seq_16S, "/home/junyuchen/Lab/16S-Prediction/"+tempFile[:-4]+"_16S_rRNA.fasta", "fasta")
                    #print(seq_16S)


                #print
            #print(record.features)
    in_handle.close()

for root,dirs,files in os.walk(dir):
    #print(dirs)
    for tempFile in files:
        if "gff" in tempFile:
            path = os.path.join(root,tempFile)
            FindrRNA(path, tempFile)
            #print(path[:-4])
            #print(dirs)
            #print(root)
            #print(tempFile[:-4])
            
            
            
