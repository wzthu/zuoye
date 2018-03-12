from BedSort import BedSort
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=BedSort(bedInput='./minidata/atac/BedForTest')

Schedule.run()
