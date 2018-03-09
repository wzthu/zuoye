from StepBase import Configure,Schedule
from FastqDump import FastqDump
from HTSeq_sam2count import HTSeq_sam2count
# Configure.setRefDir('/home/zwei/ref')
# Configure.setGenome('hg19')


HTSeq_sam2count(samInput1='./accepted_hits.sam', gtfInput1='./genome.gtf')
Schedule.run()

print('[done]')
