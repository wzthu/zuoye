# -*- coding: utf-8 -*-
"""
Created on 9th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class SortBam(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 sortOrder = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)
        self.initIO()

        self.setParam('sortOrder', sortOrder)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', os.path.join(bamOutputDir, 'aligned.sorted.bam'), '', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')

        sortOrder = self.getParam('sortOrder')

        cmdline = [
                'java -Xmx4g -jar ../../dropseq/Drop-seq_tools-1.13/jar/lib/picard-2.10.3.jar SortSam',
                'I=%s'%(bamInput[i]), 'O=%s'%(bamOutput[i]), 'SO=%s'%(sortOrder)
        ]
        self.callCmdline(cmdline)
