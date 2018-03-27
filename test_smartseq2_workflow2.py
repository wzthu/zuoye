from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSort import BamSort
from Cufflinks import Cufflinks
from Cuffmerge import Cuffmerge
from Cuffquant import Cuffquant
from Cuffnorm import Cuffnorm
from FastQC import FastQC

Configure.setIdentity('sqchen')

#Fastq-dump
fastq_dump = FastqDump(sraInput1='../../chenshengquan/zuoye/minidata/smartseq/sra')

#FastQC
fastqc = FastQC()(fastq_dump)

#Hisat2
hisat = Hisat2(ht2Idx="../../chenshengquan/zuoye/minidata/smartseq/hg19_index/genome")(fastq_dump)

sam2bam =SamToBam(threads=16)(hisat)

bamsort = BamSort()(sam2bam)

cufflinks =Cufflinks(gtfInput='../../chenshengquan/zuoye/minidata/smartseq/genome.gtf',threads=16)(bamsort)

cuffmerge = Cuffmerge(faInput1='../../chenshengquan/zuoye/minidata/smartseq/hg19.fa',gtfInput1='../../chenshengquan/zuoye/minidata/smartseq/genome.gtf',threads=16)(cufflinks)

cuffquant = Cuffquant(threads=16)(bamsort,cuffmerge)

cuffnorm = Cuffnorm(threads=16)(cuffquant,cuffmerge)

Schedule.run()

print(fastq_dump.getMarkdownEN())
print(hisat.getMarkdownEN())
print(sam2bam.getMarkdownEN())
print(bamsort.getMarkdownEN())
print(cufflinks.getMarkdownEN())
print(cuffmerge.getMarkdownEN())
print(cuffquant.getMarkdownEN())
print(cuffnorm.getMarkdownEN())