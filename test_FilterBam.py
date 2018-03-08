from FilterBam import FilterBam
from StepBase import Configure, Schedule

import os

fb = FilterBam(bamInput = './minidata/dropseq/tmp/unalign_tagged_Cell.0.bam', bamOutputDir = './minidata/dropseq/tmp', tagReject = 'XQ')
Schedule.run()
