from StepBase import Configure,Schedule
from FastqDump import FastqDump
# Configure.setRefDir('/home/zwei/ref')
# Configure.setGenome('hg19')


# adrm = AdapterRemoval(fastqInput1='chr20_1.1.fq',fastqInput2='chr20_2.1.fq')
# rs=Bowtie2()(adrm)
FastqDump(sraInput1='../bam', fastqOutputDir='./')
Schedule.run()

print('[done]')
