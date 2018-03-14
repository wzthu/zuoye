from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSort import BamSort
from Cufflinks import Cufflinks
from Cuffmerge import Cuffmerge
from Cuffquant import Cuffquant
from Cuffnorm import Cuffnorm

Configure.setIdentity('sqchen')

#Fastq-dump
fastq_dump = FastqDump(sraInput1='./minidata/smartseq/sra')

#Hisat2
hisat = Hisat2(ht2Idx="./minidata/smartseq/hg19_index/genome")(fastq_dump)

sam2bam =SamToBam(threads=16)(hisat)

bamsort = BamSort()(sam2bam)

cuffquant = Cuffquant(gtfInput='./minidata/smartseq/genome.gtf',threads=16)(bamsort)

cuffnorm = Cuffnorm(gtfInput='./minidata/smartseq/genome.gtf',threads=16)(cuffquant)

Schedule.run()