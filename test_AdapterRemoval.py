# -*- coding: utf-8 -*-
from Bowtie2 import Bowtie2
from AdapterRemoval import AdapterRemoval
from StepBase import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

adrm = AdapterRemoval(fastqInput1='./minidata/atac/end1',fastqInput2='./minidata/atac/end2')

# adrm = AdapterRemoval(fastqInput1='./minidata/atac/end1',fastqInput2='./minidata/atac/end2',
#                       fastqOutputDir1='./fastqOutputDir1',
#                       fastqOutputDir2='./fastqOutputDir2')

Schedule.run()

