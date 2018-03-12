# -*- coding: utf-8 -*-
"""

designed for scATAC-seq

@author: Weizhang

"""

from StepBase import Step


class RmDuplicates(Step):
    def __init__(self,
                 bamInput=None,
                 bamOutputDir=None,
                 picard=None,  # your picard path
                 memory='-Xmx5g',
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('picard', picard)
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('memory', memory)

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')

        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        # set all output files
        self.setOutputDir1To1('bamOutput', bamOutputDir, None, 'bam', 'bamInput')
        self.setOutputDir1To1('METRICS', bamOutputDir, None, 'picard_METRICS.txt', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bamInput', samUpstream.getOutput('bamOutput'))

    def _multiRun(self,):
        pass
        # picard = self.getParamIO('picard')
        # bamInput = self.getInputList('bamInput')
        # bamOutput = self.getOutputList('bamOutput')
        # matrixOutput = self.getOutputList('METRICS')
        #
        # for i in range(len(bamInput)):
        #     cmdline = [
        #         'java', str(self.getParam('memory')), '-jar',
        #         picard, 'MarkDuplicates REMOVE_DUPLICATES=true',
        #         'INPUT=' + bamInput[i],
        #         'OUTPUT=' + bamOutput[i],
        #         'METRICS_FILE=' + matrixOutput[i]
        #     ]
        #     result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        picard = self.getParamIO('picard')
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        matrixOutput = self.getOutputList('METRICS')
        cmdline = [
            'java', str(self.getParam('memory')), '-jar',
            picard, 'MarkDuplicates REMOVE_DUPLICATES=true',
            'INPUT=' + bamInput[i],
            'OUTPUT=' + bamOutput[i],
            'METRICS_FILE=' + matrixOutput[i]
        ]
        result = self.callCmdline('V1', cmdline)



