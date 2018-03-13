from BamToBed import BamToBed
from StepBase import Configure,Schedule

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=BamToBed(bamInput='./minidata/atac/BamForTest',
              bedOutputDir='./bedOutputDir')

Schedule.run()

