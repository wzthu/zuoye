# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/26 20:54
@Author  : Weizhang
@FileName: FragInPeak.py

count percentage of fragments in the peak region, this
program is a part of quality control for scATAC-seq
"""

from StepBase import Step, Configure
import os.path


class FragInPeak(Step):
    def __init__(self,
                 fragInput=None,
                 peakInput=None,
                 reportOutputDir=None,
                 rScript='./FragInPeak.R',
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('fragInput', fragInput)
        self.setParamIO('peakInput', peakInput)
        self.setParamIO('reportOutputDir', reportOutputDir)
        self.setParamIO('rScript', rScript)
        self.initIO()

    def impInitIO(self):
        fragInput = self.getParamIO('fragInput')
        peakInput = self.getParamIO('peakInput')
        rScript = self.getParamIO('rScript')
        reportOutputDir = self.getParamIO('reportOutputDir')

        # set all input files
        self.setInputDirOrFile('fragInput', fragInput)
        self.setInputDirOrFile('peakInput', peakInput)
        self.setInputDirOrFile('rScript', rScript)
        # set output files
        self.setOutputDir1To1('percentage', reportOutputDir, None, '_LCForSample.txt', 'peakInput', '')

        if fragInput is not None:
            self._setInputSize(len(self.getInputList('peakInput')))

    def call(self, *args):
        fragInput = args[0]
        peakInput = args[1]

        self.setParamIO('fragInput', fragInput.getOutput('bedOutput'))
        self.setParamIO('peakInput', peakInput.getOutputList('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        fragInput = self.getInputList('fragInput')
        peakInput = self.getInputList('peakInput')
        rScript = self.getInputList('rScript')
        percentage = self.getOutputList('percentage')

        cmdline = [
            'Rscript', rScript[i],
            os.path.dirname(fragInput[i]),
            peakInput[i], percentage[i]
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
