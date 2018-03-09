# -*- coding: utf-8 -*-

from BamSort import BamSort
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=BamSort(bamInput='./minidata/atac/BamForTest')

Schedule.run()

