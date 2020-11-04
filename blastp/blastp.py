from Bio.Blast.Applications import NcbiblastpCommandline

qdata = "/mnt/d/Lab/Bile-Acid/result/prodigal/GCF_000154445.1_ASM15444v1_genomic.faa"

blastpDB = "/mnt/d/Lab/Bile-Acid/db/BSH_anno.fasta"

blastpOut = "/mnt/d/Lab/Bile-Acid/result/blastp/test.xml"

blastp_cline = NcbiblastpCommandline(query=qdata, db=blastpDB, evalue=0.001,outfmt=5, out=blastpOut)

stdout, stderr = blastp_cline()