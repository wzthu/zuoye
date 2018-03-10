# -*- coding: utf-8 -*-
"""
Created on 8th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class BamToFastq(Step):
    def __init__(self,
                 bamInput = None,
                 fastqOutputDir = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('fastqOutputDir', fastqOutputDir)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        fastqOutputDir = self.getParamIO('fastqOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('fastqOutput', os.path.join(fastqOutputDir, 'unaligned_mc_tagged_polyA_filtered.fastq'), '', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        fastqOutput = self.getOutputList('fastqOutput')

        cmdline = [
                'java -Xmx4g -jar ../../dropseq/Drop-seq_tools-1.13/jar/lib/picard-2.10.3.jar SamToFastq',
                'INPUT=%s'%(bamInput[i]), 'FASTQ=%s'%(fastqOutput[i])
        ]
        self.callCmdline(cmdline)
