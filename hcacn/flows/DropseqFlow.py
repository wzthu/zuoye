# -*- coding: utf-8 -*-

from ..core import Flow, Report, Configure, Schedule
from ..steps import FastqToBam, BamMerge, TagBarcode, FilterBam, TrimAdapter, TrimPolyA, BamToFastq, StarAlign, SortBam, MergeBamAlign, TagGene, DetectError, DigitalExpression
print('222')
import os

class DropseqFlow(Flow):
    def __init__(self,
                 fastqInput1,
                 fastqInput2,
                 refDir = None,
                 refdir= None,
                 genome = None,
                 resultDir = './result',
                 cellBarcodeStart = 1,
                 cellBarcodeEnd = 16,
                 molecularBarcodeStart = 17,
                 molecularBarcodeEnd = 26,
                 adapterSeq = 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
                 starThreads = 16,
                 starOutFileNamePrefix = 'star',
                 expectedCellNum = 100):
        super(DropseqFlow, self).__init__(resultDir = resultDir,
                                          refdir = refdir,
                                          genome = genome,
                                          threads = starThreads)
        self._setParam('fastqInput1', fastqInput1)
        self._setParam('fastqInput2', fastqInput2)
        self._setParam('adapterSeq', adapterSeq)
        self._setParam('cellBarcodeStart', cellBarcodeStart)
        self._setParam('cellBarcodeEnd', cellBarcodeEnd)
        self._setParam('molecularBarcodeStart', molecularBarcodeStart)
        self._setParam('molecularBarcodeEnd', molecularBarcodeEnd)
        self._setParam('starThreads', starThreads)
        self._setParam('starOutFileNamePrefix', starOutFileNamePrefix)
        self._setParam('expectedCellNum', expectedCellNum)
        self._setParam('refDir', refDir)
        print('DropseqFlow init')

    def _call(self, *args):
        pass

    def _build(self, ):
        fastqInput1 = self._getParam('fastqInput1')
        fastqInput2 = self._getParam('fastqInput2')
        adapterSeq = self._getParam('adapterSeq')
        cellBarcodeStart = self._getParam('cellBarcodeStart')
        cellBarcodeEnd = self._getParam('cellBarcodeEnd')
        molecularBarcodeStart = self._getParam('molecularBarcodeStart')
        molecularBarcodeEnd = self._getParam('molecularBarcodeEnd')
        starThreads = self._getParam('starThreads')
        starOutFileNamePrefix = self._getParam('starOutFileNamePrefix')
        expectedCellNum = self._getParam('expectedCellNum')
        refDir = self._getParam('refDir')

        if refDir is None:
            starRef = None
            mergeRef = None
            gtfRef = None
        else:
            starRef = refDir + '/star'
            mergeRef = refDir + '/fasta'
            gtfRef = refDir + '/genes/genes.gtf'

        f2b = FastqToBam(fastqInput1 = fastqInput1, fastqInput2 = fastqInput2)
        bm = BamMerge()(f2b)
        tbc = TagBarcode(baseStart = cellBarcodeStart, baseEnd = cellBarcodeEnd,
                        barcodeRead = 1, discardRead = False,
                        tagName = 'XC', numBaseBelowQuality = 1)(bm)
        tbm = TagBarcode(baseStart = molecularBarcodeStart, baseEnd = molecularBarcodeEnd,
                        barcodeRead = 1, discardRead = True,
                        tagName = 'XM', numBaseBelowQuality = 1)(tbc)
        fb = FilterBam(tagReject = 'XQ')(tbm)
        ta = TrimAdapter(adapterSeq = adapterSeq, misMatches = 0, numBases = 5)(fb)
        tp = TrimPolyA(misMatches = 0, numBases = 6)(ta)
        b2f = BamToFastq()(tp)
        sa = StarAlign(outFileNamePrefix = starOutFileNamePrefix, genomeDir = starRef, threads = 16)(b2f)
        sb = SortBam(sortOrder = 'queryname')(sa) 
        mba = MergeBamAlign(refInputDir = mergeRef, secondAlign = False, pairedRun = False)(tp, sb)
        tg = TagGene(gtfInput = gtfRef, tag='GE')(mba)
        de = DetectError(numCells = expectedCellNum, primerSeqence = adapterSeq)(tg)
        dge = DigitalExpression(numCells = expectedCellNum)(de)
        
        rp = Report()
        rp.add('Section for FastqInput', [f2b])
        rp.add('Section for Tag Barcode', [tbc, tbm])
        rp.add('Section for Trim Reads', [ta, tp])
        rp.add('Section for STAR Alignment', [sa])
        rp.add('Section for DetectError', [de])
        rp.add('Section for DigitalExpression', [dge])

        self._setObj('FastqToBam', f2b)
        self._setObj('BamMerge', bm)
        self._setObj('TagCellBarcode', tbc)
        self._setObj('TagMolecularBarcode', tbm)
        self._setObj('FilterBam', fb)
        self._setObj('TrimAdapter', ta)
        self._setObj('TrimPolyA', tp)
        self._setObj('BamToFastq', b2f)
        self._setObj('StarAlign', sa)
        self._setObj('SortBam', sb)
        self._setObj('MergeBamAlign', mba)
        self._setObj('TagGene', tg)
        self._setObj('DetectError', de)
        self._setObj('DigitalExpression', dge)
        self._setObj('Report', rp)

        print('DropseqFlow builted')

    def _copy(self, ):
        self._linkRecursive(self._getObj('StarAlign').getOutput('logFinalOutput'),
                            self.getFinalRsDir())
        self._linkRecursive(self._getObj('StarAlign').getOutput('bamOutput'),
                            self.getFinalRsDir())
        self._linkRecursive(self._getObj('DigitalExpression').getOutput('dgeOutput'),
                            self.getFinalRsDir())
        self._linkRecursive(self._getObj('DigitalExpression').getOutput('sumOutput'),
                            self.getFinalRsDir())

