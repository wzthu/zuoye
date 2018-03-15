# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 10:58
@Author  : Weizhang
@FileName: GenPeakWithFilter.py

generate peak from macs2 summit file
summitInput must be a summit file

if topPeak > the number of peaks, R will cause a
"subscript contains out-of-bounds indices" error,
docker report "Execution halted"
"""

from StepBase import Step, Configure


class GenPeakWithFilter(Step):
    def __init__(self,
                 summitInput=None,
                 blacklist=None,
                 bedOutputDir=None,
                 overlapRate=0.2,
                 extendRange=250,
                 topPeak=50000,
                 rScript='./PeakFilter.R',
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('summitInput', summitInput)
        self.setParamIO('bedOutputDir', bedOutputDir)
        self.setParamIO('blacklist', blacklist)
        self.setParamIO('rScript', rScript)

        self.initIO()

        # set other parameters
        self.setParam('overlapRate', overlapRate)
        self.setParam('extendRange', extendRange)
        self.setParam('topPeak', topPeak)

    def impInitIO(self):
        summitInput = self.getParamIO('summitInput')
        blacklist = self.getParamIO('blacklist')
        rScript = self.getParamIO('rScript')
        bedOutputDir = self.getParamIO('bedOutputDir')

        # set all input files
        self.setInputDirOrFile('summitInput', summitInput)
        self.setInputDirOrFile('blacklist', blacklist)
        self.setInputDirOrFile('rScript', rScript)
        # set all output files
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, '_filterd.bed', 'summitInput', '')

        if summitInput is not None:
            self._setInputSize(len(self.getInputList('summitInput')))

    def call(self, *args):
        summitUpstream = args[0]

        self.setParamIO('summitInput', summitUpstream.getOutput('outputSummit'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        summitInput = self.getInputList('summitInput')
        blacklist = self.getInputList('blacklist')
        rScript = self.getInputList('rScript')
        bedOutput = self.getOutputList('bedOutput')

        cmdline = [
            'Rscript', rScript[i], summitInput[i],
            blacklist[i], bedOutput[i],
            str(self.getParam('overlapRate')), str(self.getParam('extendRange')),
            str(self.getParam('topPeak'))
        ]

        result = self.callCmdline('V1', cmdline)


