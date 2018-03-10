# -*- coding: utf-8 -*-
from SRAToFastq import SRAToFastq
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=SRAToFastq(sraInput='./minidata/atac/SraForTest')


Schedule.run()

