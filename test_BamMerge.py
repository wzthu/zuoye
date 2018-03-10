from BamMerge import BamMerge
from StepBase import Configure,Schedule

import os

bm = BamMerge(bamInput='./minidata/dropseq/tmp/', bamOutputDir='./minidata/dropseq/tmp/')
Schedule.run()
