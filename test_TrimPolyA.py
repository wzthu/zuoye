from TrimPolyA import TrimPolyA
from StepBase import Configure, Schedule

import os

ta = TrimPolyA(bamInput = './minidata/dropseq/tmp/unaligned_tagged_trimmed_smart.bam', bamOutputDir = './minidata/dropseq/tmp/',
                 sumOutputDir = './minidata/dropseq/tmp/', misMatches = 0, numBases = 5)

Schedule.run()
