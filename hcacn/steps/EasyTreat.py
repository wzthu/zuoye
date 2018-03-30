# -*- coding: utf-8 -*-

from ..core import Step, Configure
import subprocess
import os

class EasyTreat(Step):
    def __init__(self,
                 dgeInput = None,
                 outFileDir = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('dgeInput', dgeInput)
        self.setParamIO('outFileDir', outFileDir)
        #dgeOutput = None
        mtxOutput = None
        hisOutput = None
        fitOutput = None
        pkOutput = None
        if outFileDir is not None:
            #dgeOutput = os.path.join(outFileDir, 'raw.dge.txt')
            mtxOutput = os.path.join(outFileDir, 'log.tpm.with.sig.genes.txt')
            hisOutput = os.path.join(outFileDir, 'genes.detected.per.cell.eps')
            fitOutput = os.path.join(outFileDir, 'genes.cv.mu.eps')
            pkOutput = os.path.join(outFileDir, 'expr.pca.kmeans.eps')
        #self.setParamIO('dgeOutput', dgeOutput)
        self.setParamIO('mtxOutput', mtxOutput)
        self.setParamIO('hisOutput', hisOutput)
        self.setParamIO('fitOutput', fitOutput)
        self.setParamIO('pkOutput', pkOutput)
        self.initIO()

    def impInitIO(self,):
        dgeInput = self.getParamIO('dgeInput')
        outFileDir = self.getParamIO('outFileDir')
        #dgeOutput = self.getParamIO('dgeOutput')
        mtxOutput = self.getParamIO('mtxOutput')
        hisOutput = self.getParamIO('hisOutput')
        fitOutput = self.getParamIO('fitOutput')
        pkOutput = self.getParamIO('pkOutput')

        self.setInputDirOrFile('dgeInput', dgeInput)
        #self.setOutputDirNTo1('dgeOutput', dgeOutput, 'ram.dge.txt', 'dgeInput')
        self.setOutputDirNTo1('mtxOutput', mtxOutput, 'log.tpm.with.sig.genes.txt', 'dgeInput')
        self.setOutputDirNTo1('hisOutput', hisOutput, 'genes.detected.per.cell.eps', 'dgeInput')
        self.setOutputDirNTo1('fitOutput', fitOutput, 'genes.cv.mu.eps', 'dgeInput')
        self.setOutputDirNTo1('pkOutput', pkOutput, 'expr.pca.kmeans.eps', 'dgeInput')


        if dgeInput is not None:
            self._setInputSize(len(self.getInputList('dgeInput')))
        if outFileDir is None:
            self.setParamIO('outFileDir', Configure.getTmpDir())
        #if dgeOutput is None:
        #    self.setParamIO('dgeOutput', Configure.getTmpPath('ram.dge.txt'))
        if mtxOutput is None:
            self.setParamIO('mtxOutput', Configure.getTmpPath('log.tpm.with.sig.genes.txt'))
        if hisOutput is None:
            self.setParamIO('hisOutput', Configure.getTmpPath('genes.detected.per.cell.eps'))
        if fitOutput is None:
            self.setParamIO('fitOutput', Configure.getTmpPath('genes.cv.mu.eps'))
        if pkOutput is None:
            self.setParamIO('pkOutput', Configure.getTmpPath('expr.pca.kmeans.eps'))

    def call(self, *args):
        dgeUpstream = args[0]
        self.setParamIO('dgeInput', dgeUpstream.getOutput('dgeOutput'))

    def _singleRun(self, i):
        dgeInput = self.getInputList('dgeInput')
        #dgeOutput = self.getOutputList('dgeOutput')
        outFileDir = self.getParamIO('outFileDir')

        #cmdline = ['gunzip -c %s'%(dgeInput[i]), '>', dgeOutput[i]]
        #self.callCmdline('V1', cmdline)

        cmdline = ['Rscript /data/EasyTreat.R %s %s'%(dgeInput[i], outFileDir)]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## EasyTreat Result
The distribution of total genes in each cell is as follow:
![]({hisOutput})

The plot of genes' dispersion with expression:
![]({fitOutput})

The visualization of cells in cluster:
![]({pkOutput})
        """.format(hisOutput=self.getOutput('hisOutput'), fitOutput=self.getOutput('fitOutput'), pkOutput=self.getOutput('pkOutput'))
        return mdtext
