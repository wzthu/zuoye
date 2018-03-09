# -*- coding: utf-8 -*-
"""
Created on 8th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class FilterBam(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 tagReject = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)

        self.initIO()

        self.setParam('tagReject', tagReject)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', os.path.join(bamOutputDir, 'unalign_tagged_filterd.bam'), '', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')

        tagReject = self.getParam('tagReject')

        cmdline = [
                '../../dropseq/Drop-seq_tools-1.13/FilterBAM', 'INPUT=%s'%(bamInput[i]),
                'OUTPUT=%s'%(bamOutput[i]), 'TAG_REJECT=%s'%(tagReject)
                ]
        self.callCmdline(cmdline)
        
