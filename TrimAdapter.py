# -*- coding: utf-8 -*-
"""
Created on 8th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class TrimAdapter(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 sumOutputDir = None,
                 adapterSeq = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                 misMatches = 0,
                 numBases = 5,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)
        self.setParamIO('sumOutputDir', sumOutputDir)

        self.initIO()

        self.setParam('adapterSeq', adapterSeq)
        self.setParam('misMatches', misMatches)
        self.setParam('numBases', numBases)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')
        sumOutputDir = self.getParamIO('sumOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDir1To1('bamOutput', bamOutputDir, 'unaligned_tagged_trimmed_smart', 'bam', 'bamInput')
        self.setOutputDir1To1('sumOutput', sumOutputDir, 'adapter_trimming_report', 'txt', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        sumOutput = self.getOutputList('bamOutput')

        adapterSeq = self.getParam('adapterSeq')
        misMatches = self.getParam('misMatches')
        numBases = self.getParam('numBases')

        cmdline = [
                '../../dropseq/Drop-seq_tools-1.13/TrimStartingSequence',
                'INPUT=%s'%(bamInput[i]), 'OUTPUT=%s'%(bamOutput[i]),
                'OUTPUT_SUMMARY=%s'%(sumOutput[i]), 'SEQUENCE=%s'%(adapterSeq),
                'MISMATCHES=%d'%(misMatches), 'NUM_BASES=%d'%(numBases)
        ]
        self.callCmdline(cmdline)
