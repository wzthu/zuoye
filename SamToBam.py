# -*- coding: utf-8 -*-
"""

@author: Weizhang

"""

from StepBase import Step


class SamToBam(Step):
    # set default function: convert SAM to BAM
    def __init__(self,
                 samInput=None,  # <in.bam>|<in.sam>|<in.cram>
                 bamOutputDir=None,  # -o FILE  output file name [stdout]
                 threads=1,  # -@ INT number of BAM compression threads [0]
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('samInput', samInput)
        self.setParamIO('bamOutputDir', bamOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('threads', threads)

    def impInitIO(self):
        samInput = self.getParamIO('samInput')
        bamOutputDir = self.getParamIO('bamOutputDir')

        # set all input files
        self.setInputDirOrFile('samInput', samInput)
        # set all output files
        self.setOutputDir1To1('BamOutput', bamOutputDir, None, 'bam', 'samInput')

        if samInput is not None:
            self._setInputSize(len(self.getInputList('samInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('samInput', samUpstream.getOutput('samOutput'))

    def _multiRun(self,):
        samInput = self.getInputList('samInput')
        bamOutput = self.getOutputList('BamOutput')
        for i in range(len(samInput)):
            cmdline = [
                'samtools view -b -S',
                '-@', str(self.getParam('threads')),
                '-o', bamOutput[i], samInput[i]
            ]
            result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        samInput = self.getInputList('samInput')
        bamOutput = self.getOutputList('BamOutput')

        cmdline = [
            'samtools view -b -S',
            '-@', str(self.getParam('threads')),
            '-o', bamOutput[i], samInput[i]
        ]

        result = self.callCmdline(cmdline)



