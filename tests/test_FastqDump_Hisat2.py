from StepBase import Configure,Schedule
from FastqDump import FastqDump
from Hisat2 import Hisat2
# Configure.setRefDir('/home/zwei/ref')
# Configure.setGenome('hg19')


# adrm = AdapterRemoval(fastqInput1='chr20_1.1.fq',fastqInput2='chr20_2.1.fq')
# rs=Bowtie2()(adrm)
fastq_dump = FastqDump(sraInput1='../bam', fastqOutputDir='./')

fastq_dump.outputs
hisat = Hisat2(ht2Idx="../hg19_index/genome",
               samOutputDir="../smartseq2_hisat_result")(fastq_dump)

hisat.outputs
Schedule.run()

print('[done]')
