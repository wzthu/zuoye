# -*- coding: utf-8 -*-
from SRAToFastq import SRAToFastq
from StepBase import Configure,Schedule

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

# test=SRAToFastq(sraInput='./minidata/atac/SraForTest')

test=SRAToFastq(sraInput='./minidata/atac/SraForTest', fastqOutputDir='./SRAToFastq_output')

Schedule.run()

