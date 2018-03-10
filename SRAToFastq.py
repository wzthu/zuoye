# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 8:49
@Author  : Weizhang
@FileName: SRAToFastq1.py
"""


from StepBase import Step,Configure
import os


class SRAToFastq(Step):
    def __init__(self,
                 sraInput=None,
                 fastqOutputDir=None,
                 fastqc=True,  # fo fastqc or not
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('sraInput', sraInput)
        self.setParamIO('fastqOutputDir', fastqOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('fastqc', fastqc)

    def impInitIO(self):
        sraInput = self.getParamIO('sraInput')
        fastqOutputDir = self.getParamIO('fastqOutputDir')

        # set all input files
        self.setInputDirOrFile('sraInput', sraInput)
        # set all output files
        self.setOutputDir1To1('fastqOutput1', fastqOutputDir, None, '_1.fastq', 'sraInput', '')
        self.setOutputDir1To1('fastqOutput2', fastqOutputDir, None, '_2.fastq', 'sraInput', '')

        if fastqOutputDir is None:
            self.setParamIO('fastqOutputDir', Configure.getTmpDir())


        if sraInput is not None:
            self._setInputSize(len(self.getInputList('sraInput')))

    def call(self, *args):

        print("SRAToFastq has no upstream!!!\n")

    def _multiRun(self,):
        sraInput = self.getInputList('sraInput')
        fastqOutputDir = self.getParamIO('fastqOutputDir')
        fastqOutput1 = self.getOutputList('fastqOutput1')
        fastqOutput2 = self.getOutputList('fastqOutput2')

        if self.getParam('fastqc'):  # do fastqc
            for i in range(len(sraInput)):
                cmdline1 = ['fastq-dump', '--split-3',
                           sraInput[i],
                           '-O', fastqOutputDir
                           ]
                result = self.callCmdline(cmdline1)

                cmdline2 = ['fastqc', fastqOutput1[i], fastqOutput2[i]]
                result = self.callCmdline(cmdline2)

        else:  # do not do fastqc
            for i in range(len(sraInput)):
                cmdline1 = ['fastq-dump', '--split-3',
                            sraInput[i],
                            '-O', fastqOutputDir
                            ]
                result = self.callCmdline(cmdline1)



    def _singleRun(self, i):
        sraInput = self.getInputList('sraInput')
        fastqOutputDir = self.getParamIO('fastqOutputDir')
        fastqOutput1 = self.getOutputList('fastqOutput1')
        fastqOutput2 = self.getOutputList('fastqOutput2')

        if self.getParam('fastqc'):  # do fastqc
            cmdline1 = ['fastq-dump', '--split-3',
                        sraInput[i],
                        '-O', fastqOutputDir
                        ]
            result = self.callCmdline(cmdline1)

            cmdline2 = ['fastqc', fastqOutput1[i], fastqOutput2[i]]
            result = self.callCmdline(cmdline2)
        else:  # do not do fastqc
            cmdline1 = ['fastq-dump', '--split-3',
                        sraInput[i],
                        '-O', fastqOutputDir
                        ]
            result = self.callCmdline(cmdline1)


