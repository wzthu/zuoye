from DetectError import DetectError
from StepBase import Configure,Schedule

import os

de = DetectError(bamInput='./minidata/dropseq/tmp/star_gene_exon_tagged.bam', bamOutputDir='./minidata/dropseq/tmp/',
                 statsOutputDir='./minidata/dropseq/tmp/', sumOutputDir='./minidata/dropseq/tmp',
                 numCells=100, primerSeqence='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
Schedule.run()
