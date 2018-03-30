from StarAlign import StarAlign
from StepBase import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

sa = StarAlign(fastqInput='./step_00_BamToFastq/unaligned_mc_tagged_polyA_filtered.fastq', outFileNamePrefix='star',
                genomeDir = '../ref/refdata-cellranger-hg19_and_mm10-2.1.0/star/', threads=16)
Schedule.run()
