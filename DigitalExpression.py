# -*- coding: utf-8 -*-

from StepBase import Step, Configure
import subprocess
import os

class DigitalExpression(Step):
    def __init__(self,
                 bamInput = None,
                 dgeOutputDir = None,
                 sumOutputDir = None,
                 numCells = 1000,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('dgeOutputDir', dgeOutputDir)
        self.setParamIO('sumOutputDir', sumOutputDir)

        self.initIO()

        self.setParam('numCells', numCells)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        dgeOutputDir = self.getParamIO('dgeOutputDir')
        sumOutputDir = self.getParamIO('sumOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('dgeOutput', os.path.join(dgeOutputDir, 'out_gene_exon_tagged.dge.txt.gz'), '', 'bamInput')
        self.setOutputDirNTo1('sumOutput', os.path.join(sumOutputDir, 'out_gene_exon_tagged.dge.summary.txt'), '', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        dgeOutput = self.getOutputList('dgeOutput')
        sumOutput = self.getOutputList('sumOutput')

        numCells = self.getParam('numCells')

        cmdline = [
                '../../dropseq/Drop-seq_tools-1.13/DigitalExpression',
                'I=%s'%(bamInput[i]), 'O=%s'%(dgeOutput[i]), 'SUMMARY=%s'%(sumOutput[i]),
                'NUM_CORE_BARCODES=%d'%(numCells)
        ]
        self.callCmdline(cmdline)
