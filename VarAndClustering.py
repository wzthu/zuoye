# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 17:13
@Author  : Weizhang
@FileName: VarAndClustering.py
"""


from StepBase import Step, Configure
from GenPeakWithFilter import GenPeakWithFilter


class VarAndClustering(Step):
    def __init__(self,
                 bamInput=None,
                 figureOutput=None,
                 peakInput=None,
                 genome=None,
                 threads=1,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('figureOutput', figureOutput)
        self.setParamIO('peakInput', peakInput)

        self.initIO()

        # set other parameters
        self.setParam('threads', threads)
        if genome is None:
            self.setParam('genome', Configure.getGenome())
        else:
            self.setParam('genome', genome)

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        figureOutput = self.getParamIO('figureOutput')
        peakInput = self.getParamIO('peakInput')

        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        self.setInputDirOrFile('peakInput', peakInput)
        # set all output files
        self.setOutputDir1To1('var.tiff', figureOutput, None, '_var.tiff', 'peakInput', '')
        self.setOutputDir1To1('clustring.tiff', figureOutput, None, '_clustring.tiff', 'peakInput', '')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('peakInput')))

    def call(self, *args):
        bamUpstream = args[0]
        peakUpstream = args[1]
        # samOutput is from the former step (Mapping)
        self.setParamIO('bamInput', bamUpstream.getParamIO('bamOutputDir'))
        self.setParamIO('peakInput', peakUpstream.getOutputList('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bamInput = self.getParamIO('bamInput')
        peakInput = self.getParamIO('peakInput')
        varOutput = self.getOutputList('var.tiff')
        clustringOutput = self.getOutputList('clustring.tiff')
        cmdline = [
            'Rscript ChromVarUtility.R',
            bamInput,
            peakInput,
            str(self.getParam('threads')),
            str(self.getParam('genome')),
            varOutput[i], clustringOutput[i]
        ]
        result = self.callCmdline(cmdline)

