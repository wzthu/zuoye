# -*- coding: utf-8 -*-
"""

designed for scATAC-seq

@author: Weizhang

call peaks using macs2

"""

from StepBase import Step, Configure
from pathlib import Path


class PeakCalling(Step):
    def __init__(self,
                 bedInput=None,
                 outputDir=None,
                 format='BED',
                 genome='hs',
                 nomodel=True,
                 nolambda=True,
                 keepdup='--keep-dup all',  # can be set to ''
                 callsummits=True,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bedInput', bedInput)
        self.setParamIO('outputDir', outputDir)

        self.initIO()

        self.setParam('format', format)
        self.setParam('genome', genome)
        self.setParam('nomodel', nomodel)
        self.setParam('nolambda', nolambda)
        self.setParam('keepdup', keepdup)
        self.setParam('callsummits', callsummits)

    def impInitIO(self):
        bedInput = self.getParamIO('bedInput')
        outputDir = self.getParamIO('outputDir')

        # set all input files
        self.setInputDirOrFile('bedInput', bedInput)
        # set all output files
        self.setOutputDir1To1('outputSummit', outputDir, None, '_summits.bed', 'bedInput', '')
        self.setOutputDir1To1('narrowPeak', outputDir, None, '_peaks.narrowPeak', 'bedInput', '')
        self.setOutputDir1To1('outputXls', outputDir, None, '_peaks.xls', 'bedInput', '')


        if bedInput is not None:
            self._setInputSize(len(self.getInputList('bedInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bedInput', samUpstream.getOutput('bedOutput'))

    def _multiRun(self,):
        pass
        # bedInput = self.getInputList('bedInput')
        # bedOutput = self.getOutputList('bedOutput')
        # for i in range(len(bedInput)):
        #     cmdline = [
        #         'sortBed -i',
        #         bedInput[i],
        #         '>',
        #         bedOutput[i]
        #     ]
        #     result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        bedInput = self.getInputList('bedInput')
        outputSummit = self.getOutputList('outputSummit')
        narrowPeak = self.getOutputList('narrowPeak')
        outputXls = self.getOutputList('outputXls')
        cmdline = [
            'macs2 callpeak -t', bedInput[i],
            '-f', str(self.getParam('format')),
            '-g', str(self.getParam('genome')),
            '-n', Path(bedInput[i]).stem,
            self.getBoolParamCmd('--nomodel', 'nomodel'),
            self.getBoolParamCmd('--nomodel', 'nomodel'),
            self.getBoolParamCmd('--nolambda', 'nolambda'),
            self.getParam('keepdup'),
            self.getBoolParamCmd('--call-summits', 'callsummits'),
            '&&',
            'mv', './MergedAllFiles_summits.bed', outputSummit[i],
            '&&',
            'mv', './MergedAllFiles_peaks.narrowPeak', narrowPeak[i],
            '&&',
            'mv', './MergedAllFiles_peaks.xls', outputXls[i]
        ]
        result = self.callCmdline('V2', cmdline)

