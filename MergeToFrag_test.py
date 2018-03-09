# -*- coding: utf-8 -*-
from MergeToFrag import MergeToFrag
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

test=MergeToFrag(bedInput='./minidata/atac/BedForTest')

Schedule.run()