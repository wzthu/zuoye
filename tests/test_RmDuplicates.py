# -*- coding: utf-8 -*-
from RmDuplicates import RmDuplicates
from StepBase import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=RmDuplicates(bamInput='./minidata/atac/BamForTest',
                  memory='-Xmx4g')

Schedule.run()

