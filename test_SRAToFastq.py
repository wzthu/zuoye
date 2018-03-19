# -*- coding: utf-8 -*-
from SRAToFastq import SRAToFastq
from StepBase import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=SRAToFastq(sraInput='./minidata/atac/SraForTest')

Schedule.run()

f = open("/data8t_1/zhangwei1/11111111.Rmd", "w")
f.write(test.getMarkdownEN())
f.close()


