from FastqToBam import FastqToBam
from StepBase import Configure, Schedule

lcy_b2s = FastqToBam(fastqInput1='./minidata/dropseq/read1', fastqInput2='./minidata/dropseq/read2', bamOutputDir='./minidata/dropseq/tmp1')
Schedule.run()

print('')
