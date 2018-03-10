from DigitalExpression import DigitalExpression
from StepBase import Configure,Schedule

import os

de = DigitalExpression(bamInput='./minidata/dropseq/tmp/out_gene_exon_tagged.bam', dgeOutputDir='./minidata/dropseq/tmp/',
                 sumOutputDir='./minidata/dropseq/tmp', numCells=100)
Schedule.run()
