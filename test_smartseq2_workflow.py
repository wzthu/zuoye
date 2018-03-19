from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSort import BamSort
from Cufflinks import Cufflinks
from Cuffmerge import Cuffmerge
from Cuffquant import Cuffquant
from Cuffdiff import Cuffdiff

Configure.setIdentity('sqchen')

#Fastq-dump
fastq_dump = FastqDump(sraInput1='../../chenshengquan/zuoye/minidata/smartseq/sra')

#Hisat2
hisat = Hisat2(ht2Idx="../../chenshengquan/zuoye/minidata/smartseq/hg19_index/genome")(fastq_dump)

# Bam2Sam
sam2bam =SamToBam(threads=16)(hisat)

# #BamSort
bamsort = BamSort()(sam2bam)

# #Cufflinks
cufflinks =Cufflinks(gtfInput='../../chenshengquan/zuoye/minidata/smartseq/genome.gtf',threads=16)(bamsort)

cuffmerge=Cuffmerge(faInput1='../../chenshengquan/zuoye/minidata/smartseq/hg19.fa',gtfInput1='../../chenshengquan/zuoye/minidata/smartseq/genome.gtf',threads=16)(cufflinks)

cuffquant = Cuffquant(threads=16)(bamsort,cuffmerge)

cuffdiff = Cuffdiff(faInput='../../chenshengquan/zuoye/minidata/smartseq/hg19.fa',threads=16)(cuffmerge,cuffquant)

Schedule.run()

print(fastq_dump.getMarkdownEN())
print(hisat.getMarkdownEN())
print(sam2bam.getMarkdownEN())
print(bamsort.getMarkdownEN())
print(cufflinks.getMarkdownEN())
print(cuffmerge.getMarkdownEN())
print(cuffquant.getMarkdownEN())
print(cuffdiff.getMarkdownEN())