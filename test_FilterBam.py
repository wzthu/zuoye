from FilterBam import FilterBam
from StepBase import Configure, Schedule

import os

fb = BamMerge(bamInput = './minidata/dropseq/tmp/unaligned_tagged_CellMolecular.bam', bamOutputDir = './minidata/dropseq/tmp')
Schedule.run()
