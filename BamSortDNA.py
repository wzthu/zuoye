# -*- coding: utf-8 -*-
"""

designed for scATAC-seq

@author: Weizhang

"""

from StepBase import Step


class BamSort(Step):
    def __init__(self,
                 bamInput=None,
                 bamOutputDir=None,
                 threads=1,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('threads', threads)

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')

        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        # set all output files
        self.setOutputDir1To1('BamOutput', bamOutputDir, None, 'bam', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bamInput', samUpstream.getOutput('BamOutput'))

    def _multiRun(self,):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('BamOutput')

        for i in range(len(bamInput)):
            cmdline = [
                'samtools sort -O BAM',
                '-@', str(self.getParam('threads')),
                '-o', bamOutput[i], bamInput[i]
            ]
            result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        samInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('BamOutput')

        cmdline = [
            'samtools sort -O BAM',
            '-@', str(self.getParam('threads')),
            '-o', bamOutput[i], samInput[i]
        ]

        result = self.callCmdline(cmdline)



