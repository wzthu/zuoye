from BedSort import BedSort
from StepBase import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=BedSort(bedInput='./minidata/atac/BedForTest')

Schedule.run()
