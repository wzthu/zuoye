# -*- coding: utf-8 -*-
from RmDuplicates import RmDuplicates
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=RmDuplicates(bamInput='./minidata/atac/BamForTest',
                  picard='/home/wzhang/software/Picard/picard.jar',
                  memory='-Xmx4g',
                  bamOutputDir='/home/wzhang/scATAC_test/pipeline_test/111')

Schedule.run()