from SortBam import SortBam
from StepBase import Configure,Schedule

import os

sb = SortBam(bamInput='./minidata/dropseq/tmp/starAligned.out.sam', bamOutputDir='./minidata/dropseq/tmp/', sortOrder = 'queryname')
Schedule.run()
