# -*- coding: utf-8 -*-
"""
Created on 9th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class TagGene(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutputDir = None,
                 gtfInput = None,
                 tag = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)
        self.setParamIO('gtfInput', gtfInput)

        self.initIO()

        self.setParam('tag', tag)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')
        gtfInput = self.getParamIO('gtfInput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setInputDirOrFile('gtfInput', gtfInput)
        self.setOutputDirNTo1('bamOutput', os.path.join(bamOutputDir, 'star_gene_exon_tagged.bam'), '', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamInput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        gtfInput = self.getInputList('gtfInput')
        bamOutput = self.getOutputList('bamOutput')

        tag = self.getParam('tag')

        cmdline = [
                '../../dropseq/Drop-seq_tools-1.13/TagReadWithGeneExon', 'I=%s'%(bamInput[i]),
                'O=%s'%(bamOutput[i]), 'ANNOTATIONS_FILE=%s'%(gtfInput[i]), 'TAG=%s'%(tag)
        ]
        self.callCmdline(cmdline)
