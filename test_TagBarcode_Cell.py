from TagBarcode import TagBarcode
from StepBase import Configure,Schedule

import os

bm = TagBarcode(bamInput='./minidata/dropseq/tmp/unaligned_data.bam', bamOutputDir='./minidata/dropseq/tmp/',
                sumOutputDir='./minidata/dropseq/tmp/', baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)
Schedule.run()
