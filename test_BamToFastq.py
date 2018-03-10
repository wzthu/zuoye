from BamToFastq import BamToFastq
from StepBase import Configure, Schedule

import os

b2f = BamToFastq(bamInput = './minidata/dropseq/tmp/unaligned_mc_tagged_polyA_filtered.bam', fastqOutputDir = './minidata/dropseq/tmp/')

Schedule.run()
