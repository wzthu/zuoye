# -*- coding: utf-8 -*-
"""
Created on Six MAR

@author: CyLiu
"""

from StepBase import Step,Configure
import subprocess
import os

class FastqToBam(Step):
    def __init__(self,
                 fastqInput1 = None,
                 fastqInput2 = None,
                 bamOutputDir = None,
                 #mapRsOutputDir = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('fastqInput1', fastqInput1)
        self.setParamIO('fastqInput2', fastqInput2)
        self.setParamIO('bamOutputDir', bamOutputDir)
        #self.setParamIO('mapRsOutputDir', mapRsOutputDir)

        self.initIO()


    def impInitIO(self,):
        fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput2 = self.getParamIO('fastqInput2')
        bamOutputDir = self.getParamIO('bamOutputDir')
        #mapRsOutputDir = self.getParamIO('mapRsOutputDir')

        self.setInputDirOrFile('fastqInput1', fastqInput1)
        self.setInputDirOrFile('fastqInput2', fastqInput2)

        self.setOutputDir1To1('bamOutput', bamOutputDir, 'sample', 'bam', 'fastqInput1')
        #self.setOutputDir1To1('mapRsOutput', mapRsOutputDir, 'sample', 'txt', 'fastqInput1')
        if fastqInput1 is not None:
            self._setInputSize(len(self.getInputList('fastqInput1')))

    def call(self, *args):
        fastqUpstream = args[0]
        self.setParamIO('fastqInput1', fastqUpstream.getOutput('fastqOutput1'))
        self.setParamIO('fastqInput2', fastqUpstream.getOutput('fastqOutput2'))

    def _singleRun(self, i):
        fastqInput1 = self.getInputList('fastqInput1')
        fastqInput2 = self.getInputList('fastqInput2')
        bamOutput = self.getOutputList('bamOutput')
        #mapRsOutput = self.getOutputList('mapRsOutput')
        cmdline = ['java -jar ../../dropseq/Drop-seq_tools-1.13/jar/lib/picard-2.10.3.jar FastqToSam',
                  'F1=%s F2=%s O=%s SM=for_tool_testing'%(fastqInput1[i], fastqInput2[i], bamOutput[i])]
        self.callCmdline(cmdline)
