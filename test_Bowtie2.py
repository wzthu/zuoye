from Bowtie2 import Bowtie2
from AdapterRemoval import AdapterRemoval
from StepBase import Configure,Schedule

import os

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

rs=Bowtie2(fastqInput1='./minidata/atac/end1',
           fastqInput2='./minidata/atac/end2',
           samOutputDir='./samOutputDir')

Schedule.run()

#Schedule.stopDocker('V1')


