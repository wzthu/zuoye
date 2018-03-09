from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSort import BamSort
from Cufflinks import Cufflinks


#Fastq-dump
fastq_dump = FastqDump(sraInput1='../bam', fastqOutputDir='./')
print(fastq_dump.outputs)

#Hisat2
hisat = Hisat2(ht2Idx="../hg19_index/genome",
               samOutputDir="../smartseq2_hisat_result")(fastq_dump)
print (hisat.outputs)

# Bam2Sam
Configure.setRefDir('../hg19_bowtie2/')
Configure.setGenome('hg19')
sam2bam =SamToBam(threads=5)(hisat)
print (sam2bam.outputs)

#BamSort
bamsort = BamSort()(sam2bam)
print (bamsort.outputs)

#Cufflinks

cufflinks =Cufflinks(gtfInput=['../genome.gtf'],outputDir=['./'])(bamsort)
print (cufflinks.outputs)

Schedule.run()