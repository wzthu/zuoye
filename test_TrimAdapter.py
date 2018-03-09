from TrimAdapter import TrimAdapter
from StepBase import Configure, Schedule

import os

ta = TrimAdapter(bamInput = './minidata/dropseq/tmp/unalign_tagged_filterd.bam', bamOutputDir = './minidata/dropseq/tmp/',
                 sumOutputDir = './minidata/dropseq/tmp/', adapterSeq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', misMatches = 0, numBases = 5)

Schedule.run()
