# -*- coding: utf-8 -*-
from StepBase import Step, Configure
import subprocess
import os

class DetectError(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 statsOutputDir = None,
                 sumOutputDir = None,
                 numCells = 1000,
                 primerSeqence = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)
        self.setParamIO('statsOutputDir', statsOutputDir)
        self.setParamIO('sumOutputDir', sumOutputDir)

        self.initIO()

        self.setParam('numCells', numCells)
        self.setParam('primerSeqence', primerSeqence)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')
        statsOutputDir = self.getParamIO('statsOutputDir')
        sumOutputDir = self.getParamIO('sumOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', os.path.join(bamOutputDir, 'out_gene_exon_tagged.bam'), '', 'bamInput')
        self.setOutputDirNTo1('statsOutput', os.path.join(statsOutputDir, 'synthesis_stats.txt'), '', 'bamInput')
        self.setOutputDirNTo1('sumOutput', os.path.join(sumOutputDir, 'synthesis_stats.summary.txt'), '', 'bamInput')
        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        statsOutput = self.getOutputList('statsOutput')
        sumOutput = self.getOutputList('sumOutput')

        numCells = self.getParam('numCells')
        primerSeqence = self.getParam('primerSeqence')

        cmdline = [
                '../../dropseq/Drop-seq_tools-1.13/DetectBeadSynthesisErrors',
                'I=%s'%(bamInput[i]), 'O=%s'%(bamOutput[i]), 'OUTPUT_STATS=%s'%(statsOutput[i]),
                'SUMMARY=%s'%(sumOutput[i]), 'NUM_BARCODES=%d'%(4*numCells),
                'PRIMER_SEQUENCE=%s'%(primerSeqence)
        ]
        self.callCmdline(cmdline)
