# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 10:58
@Author  : Weizhang
@FileName: GenPeakWithFilter.py

generate peak from macs2 summit file
summitInput must be a summit file
"""

from StepBase import Step, Configure


class GenPeakWithFilter(Step):
    def __init__(self,
                 summitInput=None,
                 blacklist=None,
                 bedOutputDir=None,
                 overlapRate=0.2,
                 extendRange=250,
                 topPeak=0,
                 rScript='./PeakFilter.R',
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('summitInput', summitInput)
        self.setParamIO('bedOutputDir', bedOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('blacklist', blacklist)
        self.setParam('overlapRate', overlapRate)
        self.setParam('extendRange', extendRange)
        self.setParam('topPeak', topPeak)
        self.setParam('rScript', rScript)

    def impInitIO(self):
        summitInput = self.getParamIO('summitInput')
        bedOutputDir = self.getParamIO('bedOutputDir')

        # set all input files
        self.setInputDirOrFile('summitInput', summitInput)
        # set all output files
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, '_filterd.bed', 'summitInput', '')

        if summitInput is not None:
            self._setInputSize(len(self.getInputList('summitInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('summitInput', samUpstream.getOutput('summitOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        summitInput = self.getInputList('summitInput')
        bedOutput = self.getOutputList('bedOutput')

        cmdline = [
            'Rscript', str(self.getParam('rScript')), summitInput[i],
            str(self.getParam('blacklist')), bedOutput[i],
            str(self.getParam('overlapRate')), str(self.getParam('extendRange')),
            str(self.getParam('topPeak'))
        ]

        result = self.callCmdline(cmdline)


