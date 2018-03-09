from FastqToBam import FastqToBam
from BamMerge import BamMerge
from TagBarcode import TagBarcode
from FilterBam import FilterBam
from TrimAdapter import TrimAdapter
from TrimPolyA import TrimPolyA
from BamToFastq import BamToFastq
from StarAlign import StarAlign
from SortBam import SortBam
from MergeBamAlign import MergeBamAlign
from TagGene import TagGene
from DetectError import DetectError
from DigitalExpression import DigitalExpression
from StepBase import Configure,Schedule

import os
f2b = FastqToBam(fastqInput1 = './minidata/dropseq/read1', fastqInput2 = './minidata/dropseq/read2', bamOutputDir = './minidata/dropseq/tmp1')
bm = BamMerge(bamOutputDir = './minidata/dropseq/tmp/merged.bam')(f2b)
tbc = TagBarcode(bamOutputDir = './minidata/dropseq/tmp', sumOutputDir = './minidata/dropseq/tmp', baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)(bm)
tbm = TagBarcode(bamOutputDir = './minidata/dropseq/tmp', sumOutputDir = './minidata/dropseq/tmp', baseStart = 17, baseEnd = 26, baseQuality = 10,
                barcodeRead = 1, discardRead = True, tagName = 'XM', numBaseBelowQuality = 1)(tbc)
fb = FilterBam(bamOutputDir = './minidata/dropseq/tmp', tagReject = 'XQ')(tbm)
ta = TrimAdapter(bamOutputDir = './minidata/dropseq/tmp/', sumOutputDir = './minidata/dropseq/tmp/', adapterSeq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
                  misMatches = 0, numBases = 5)(fb)
tp = TrimPolyA(bamOutputDir = './minidata/dropseq/tmp/', sumOutputDir = './minidata/dropseq/tmp/', misMatches = 0, numBases = 6)(ta)
b2f = BamToFastq(fastqOutputDir = './minidata/dropseq/tmp/')(tp)
sa = StarAlign(outFileDir='./minidata/dropseq/tmp/star', genomeDir = '../../dropseq/refdata-cellranger-hg19_and_mm10-2.1.0/star/', threads=1)(b2f)
sb = SortBam(bamOutputDir='./minidata/dropseq/tmp/', sortOrder = 'queryname')(sa)
mba = MergeBamAlign(bamOutputDir='./minidata/dropseq/tmp', refSequence='../../dropseq/refdata-cellranger-hg19_and_mm10-2.1.0/fasta/genome.fa',
                secondAlign=False, pairedRun=False)(tp, sb)
tg = TagGene(bamOutputDir='./minidata/dropseq/tmp/', gtfInput = '../../dropseq/refdata-cellranger-hg19_and_mm10-2.1.0/genes/genes.gtf', tag='GE')(mba)
de = DetectError(bamOutputDir='./minidata/dropseq/tmp/', statsOutputDir='./minidata/dropseq/tmp/', sumOutputDir='./minidata/dropseq/tmp',
                 numCells=100, primerSeqence='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')(tg)
dge = DigitalExpression(dgeOutputDir='./minidata/dropseq/tmp/', sumOutputDir='./minidata/dropseq/tmp', numCells=100)(de)
Schedule.run()
