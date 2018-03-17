from BedSort import BedSort
from StepBase import Configure,Schedule

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=BedSort(bedInput='./minidata/atac/BedForTest',
             bedOutputDir='./bedOutputDir')

Schedule.run()
