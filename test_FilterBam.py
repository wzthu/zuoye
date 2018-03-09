from FilterBam import FilterBam
from StepBase import Configure, Schedule

import os

fb = FilterBam(bamInput = './minidata/dropseq/tmp/unalign_tagged_XM.bam', bamOutputDir = './minidata/dropseq/tmp', tagReject = 'XQ')
Schedule.run()
