from TagBarcode import TagBarcode
from StepBase import Configure,Schedule

import os

bm = TagBarcode(bamInput='./minidata/dropseq/tmp/unalign_tagged_XC.bam', bamOutputDir='./minidata/dropseq/tmp',
                sumOutputDir='./minidata/dropseq/tmp', baseStart = 17, baseEnd = 26, baseQuality = 10,
                barcodeRead = 1, discardRead = True, tagName = 'XM', numBaseBelowQuality = 1)
Schedule.run()
