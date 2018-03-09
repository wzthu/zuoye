# -*- coding: utf-8 -*-
"""
Created on 7th Mar

@author: CyLiu
"""

from StepBase import Step,Configure
import subprocess
import os

class BamMerge(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        #self.setOutputDir1To1('bamOutput', bamOutputDir, 'merged', 'bam', 'bamInput')
        self.setOutputDirNTo1('bamOutput', bamOutputDir, 'merged.bam', 'bamInput')
        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))
        self._setMultiRun()

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _multiRun(self):
        bamInput = self.getInputList('bamInput')
        bamoOutput = self.getOutputList('bamOutput')
        cmdline = ['samtools merge', bamoOutput[0]]
        for i in range(len(bamInput)):
            cmdline.append(bamInput[i])
        result = self.callCmdline(cmdline)
