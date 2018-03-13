# -*- coding: utf-8 -*-

from BamSort import BamSort
from StepBase import Configure,Schedule

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=BamSort(bamInput='./minidata/atac/BamForTest',
             bamOutputDir='./bamOutputDir')

Schedule.run()

