from StarAlign import StarAlign
from StepBase import Configure,Schedule

import os

sa = StarAlign(fastqInput='./minidata/dropseq/tmp/unaligned_mc_tagged_polyA_filtered.fastq', outFileDir='./minidata/dropseq/tmp/star',
                genomeDir = '../../dropseq/refdata-cellranger-hg19_and_mm10-2.1.0/star/', threads=1)
Schedule.run()
