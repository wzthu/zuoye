from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSortRNA import Bamsort
from Cufflinks import Cufflinks


#Fastq-dump
fastq_dump = FastqDump(sraInput1='../bam', fastqOutputDir='./')
print(fastq_dump.outputs)

#Hisat2
hisat = Hisat2(ht2Idx="../../yinqijin/hg19_index/genome",
               samOutputDir="../smartseq2_hisat_result")(fastq_dump)
print (hisat.outputs)

# Bam2Sam
Configure.setRefDir('../../yinqijin/hg19_bowtie2/')
Configure.setGenome('hg19')
sam2bam =SamToBam(threads=5)(hisat)
print (sam2bam.outputs)

#BamSort
bamsort = Bamsort()(sam2bam)
print (bamsort.outputs)

#Cufflinks

cufflinks =Cufflinks(gtfInput='../../yinqijin/genome.gtf',outputDir='../smartseq2_cufflinks_result')(bamsort)
print (cufflinks.outputs)

Schedule.run()