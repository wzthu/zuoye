from FastqToBam import FastqToBam
from BamMerge import BamMerge
from TagBarcode import TagBarcode
from StepBase import Configure,Schedule

import os
fb = FastqToBam(fastqInput1 = './minidata/dropseq/read1', fastqInput2 = './minidata/dropseq/read2', bamOutputDir = '../minidata/dropseq/tmp2')
bm = BamMerge(bamOutputDir = './minidata/dropseq/tmp3/merged.bam')(fb)
tb = TagBarcode(bamOutputDir = './minidata/dropseq/tmp4', sumOutputDir = './minidata/dropseq/tmp4', baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)(bm)
Schedule.run()
