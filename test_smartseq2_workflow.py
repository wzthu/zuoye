from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSort import BamSort
from Cufflinks import Cufflinks
from Cuffmerge import Cuffmerge

Configure.setIdentity('sqchen')

#Fastq-dump
fastq_dump = FastqDump(sraInput1='./minidata/smartseq/sra')
# print(fastq_dump.outputs)

#Hisat2
hisat = Hisat2(ht2Idx="./minidata/smartseq/hg19_index/genome")(fastq_dump)
# print (hisat.outputs)

# Bam2Sam
# Configure.setRefDir('../../yinqijin/hg19_bowtie2/')
# Configure.setGenome('hg19')
sam2bam =SamToBam(threads=5)(hisat)
# print (sam2bam.outputs)

# #BamSort
bamsort = BamSort()(sam2bam)
# print (bamsort.outputs)

# #Cufflinks
cufflinks =Cufflinks(gtfInput='./minidata/smartseq/genome.gtf',threads=16)(bamsort)
# print (cufflinks.outputs)
cuffmerge=Cuffmerge(faInput1='./minidata/smartseq/hg19.fa',gtfInput1='./minidata/smartseq/genome.gtf',threads=16)(cufflinks)
Schedule.run()