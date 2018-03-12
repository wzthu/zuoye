# -*- coding: utf-8 -*-
from SamToBam import SamToBam
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=SamToBam(samInput='./minidata/atac/SamForTest',
              threads=5)

Schedule.run()