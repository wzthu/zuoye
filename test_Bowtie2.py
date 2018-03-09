from Bowtie2 import Bowtie2
from AdapterRemoval import AdapterRemoval
from StepBase import Configure,Schedule

import os

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

rs=Bowtie2(fastqInput1='./minidata/atac/end1',
           fastqInput2='./minidata/atac/end2')

Schedule.run()


