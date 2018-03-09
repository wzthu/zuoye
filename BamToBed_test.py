from BamToBed import BamToBed
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=BamToBed(bamInput='./minidata/atac/BamForTest')

Schedule.run()

