# -*- coding: utf-8 -*-
"""
Created on 7th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class TagBarcode(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 sumOutputDir = None,
                 baseStart = 1,
                 baseEnd = 12,
                 baseQuality = 10,
                 barcodedRead = 1,
                 discardRead = False,
                 tagName = None,
                 numBaseBelowQuality = 1,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)
        self.setParamIO('sumOutputDir', sumOutputDir)
        self.initIO()

        self.setParam('baseStart', baseStart)
        self.setParam('baseEnd', baseEnd)
        self.setParam('baseQuality', baseQuality)
        self.setParam('barcodedRead', barcodedRead)
        self.setParam('discardRead', discardRead)
        self.setParam('tagName', tagName)
        self.setParam('numBaseBelowQuality', numBaseBelowQuality)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')
        sumOutputDir = self.getParamIO('sumOutputDir')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDir1To1('bamOutput', bamOutputDir, 'unalign_tagged_Cell', 'bam', 'bamInput')
        self.setOutputDir1To1('sumOutput', sumOutputDir, 'unalign_tagged_Cellular.bam_summary', 'txt', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        sumOutput = self.getOutputList('sumOutput')

        baseStart = self.getParam('baseStart')
        baseEnd = self.getParam('baseEnd')
        baseQuality = self.getParam('baseQuality')
        barcodedRead = self.getParam('barcodedRead')
        discardRead = self.getParam('discardRead')
        tagName = self.getParam('tagName')
        numBaseBelowQuality = self.getParam('numBaseBelowQuality')

        cmdline = [
                '../../dropseq/Drop-seq_tools-1.13/TagBamWithReadSequenceExtended',
                'INPUT=%s'%(bamInput[i]), 'OUTPUT=%s'%(bamOutput[i]), 'SUMMARY=%s'%(sumOutput[i]),
                'BASE_RANGE=%d-%d'%(baseStart, baseEnd), 'BASE_QUALITY=%d'%(baseQuality),
                'BARCODED_READ=%d'%(barcodedRead), 'DISCARD_READ=%s'%(str(discardRead)),
                'TAG_NAME=%s'%(tagName), 'NUM_BASES_BELOW_QUALITY=%d'%(numBaseBelowQuality)
                ]
        self.callCmdline(cmdline)
