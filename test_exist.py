from FastqToBam import FastqToBam
from BamMerge import BamMerge
from TagBarcode import TagBarcode
from FilterBam import FilterBam
from StepBase import Configure,Schedule

import os
fb = FastqToBam(fastqInput1 = './minidata/dropseq/read1', fastqInput2 = './minidata/dropseq/read2', bamOutputDir = './minidata/dropseq/tmp2')
bm = BamMerge(bamOutputDir = './minidata/dropseq/tmp3/merged.bam')(fb)
tbc = TagBarcode(bamOutputDir = './minidata/dropseq/tmp4', sumOutputDir = './minidata/dropseq/tmp4', baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)(bm)
tbm = TagBarcode(bamOutputDir = './minidata/dropseq/tmp5', sumOutputDir = './minidata/dropseq/tmp5', baseStart = 17, baseEnd = 26, baseQuality = 10,
                barcodeRead = 1, discardRead = True, tagName = 'XM', numBaseBelowQuality = 1)(tbc)
fb = FilterBam(bamOutputDir = './minidata/dropseq/tmp6', tagReject = 'XQ')(tbm)
Schedule.run()
