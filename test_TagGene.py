from TagGene import TagGene
from StepBase import Configure,Schedule

import os

tg = TagGene(bamInput='./minidata/dropseq/tmp/merged.bam', bamOutputDir='./minidata/dropseq/tmp/',
                gtfInput = '../../dropseq/refdata-cellranger-hg19_and_mm10-2.1.0/genes/genes.gtf',
                tag='GE')
Schedule.run()
