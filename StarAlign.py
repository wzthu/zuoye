# -*- coding: utf-8 -*-
"""
Created on 9th Mar

@author: CyLiu
"""

from StepBase import Step, Configure
import subprocess
import os

class StarAlign(Step):
    def __init__(self,
                 fastqInput = None,
                 outFileDir = None,
                 genomeDir = None,
                 threads = None,
                 #outSamType = 'BAM'
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, ** kwargs)

        self.setParamIO('fastqInput', fastqInput)
        self.setParamIO('outFileDir', outFileDir)
        self.setParamIO('genomeDir', genomeDir)

        self.initIO()

        #self.setParam('outSamType', outSamType)
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads', threads)

    def impInitIO(self,):
        fastqInput = self.getParamIO('fastqInput')
        outFileDir = self.getParamIO('outFileDir')
        genomeDir = self.getParamIO('genomeDir')

        for i in ['chrLength.txt', 'chrName.txt',
                  'exonGeTrInfo.tab', 'geneInfo.tab',
                  'genomeParameters.txt', 'SAindex',
                  'sjdbList.fromGTF.out.tab', 'transcriptInfo.tab',
                  'chrNameLength.txt', 'chrStart.txt',
                  'exonInfo.tab', 'Genome', 'SA',
                  'sjdbInfo.txt', 'sjdbList.out.tab']:
            self.setInputDirOrFile(i, os.path.join(genomeDir, i))

        self.setInputDirOrFile('fastqInput', fastqInput)
        self.setOutputDirNTo1('bamOutput', os.path.join(outFileDir, 'starAligned.out.bam'), '', 'fastqInput')
        
        if fastqInput is not None:
            self._setInputSize(len(self.getInputList('fastqInput')))

    def call(self, *args):
        fastqUpstream = args[0]
        self.setParamIO('fastqInput', fastqUpstream.getOutput('fastqOutput'))

    def _singleRun(self, i):
        fastqInput = self.getInputList('fastqInput')
        outFileDir = self.getParamIO('outFileDir')
        genomeDir = self.getParamIO('genomeDir')
        
        threads = self.getParam('threads')
        cmdline = [
                '../../dropseq/STAR/source/STAR', '--runThreadN %d'%(threads),'--genomeDir %s'%(genomeDir),
                '--readFilesIn %s'%(fastqInput[i]), '--outFileNamePrefix %s'%(outFileDir)
        ]
        self.callCmdline(cmdline)
