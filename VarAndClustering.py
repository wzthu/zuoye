# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 17:13
@Author  : Weizhang
@FileName: VarAndClustering.py
"""


from StepBase import Step, Configure
from GenPeakWithFilter import GenPeakWithFilter
import os.path


class VarAndClustering(Step):
    def __init__(self,
                 bamInput=None,
                 figureOutput=None,
                 peakInput=None,
                 genome=None,
                 threads=1,
                 rScript='./ChromVarUtility.R',
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('figureOutput', figureOutput)
        self.setParamIO('peakInput', peakInput)
        self.setParamIO('rScript', rScript)

        self.initIO()

        # set other parameters
        self.setParam('threads', threads)
        if genome is None:
            self.setParam('genome', Configure.getGenome())
        else:
            self.setParam('genome', genome)

        self._setUpstreamSize(2)

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        figureOutput = self.getParamIO('figureOutput')
        peakInput = self.getParamIO('peakInput')
        rScript = self.getParamIO('rScript')

        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        self.setInputDirOrFile('peakInput', peakInput)
        self.setInputDirOrFile('rScript', rScript)
        # set all output files
        self.setOutputDir1To1('var.tiff', figureOutput, None, '_var.tiff', 'peakInput', '')
        self.setOutputDir1To1('clustring.tiff', figureOutput, None, '_clustring.tiff', 'peakInput', '')

        if peakInput is not None:
            self._setInputSize(len(self.getInputList('peakInput')))

    def call(self, *args):
        bamUpstream = args[0]
        peakUpstream = args[1]

        self.setParamIO('bamInput', bamUpstream.getOutputList('bamOutput'))
        self.setParamIO('peakInput', peakUpstream.getOutputList('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        peakInput = self.getInputList('peakInput')
        rScript = self.getInputList('rScript')
        varOutput = self.getOutputList('var.tiff')
        clustringOutput = self.getOutputList('clustring.tiff')
        cmdline = [
            'Rscript', rScript[i],
            os.path.dirname(bamInput[i]),
            peakInput[i],
            str(self.getParam('threads')),
            str(self.getParam('genome')),
            varOutput[i], clustringOutput[i]
        ]
        result = self.callCmdline('V1', cmdline)

