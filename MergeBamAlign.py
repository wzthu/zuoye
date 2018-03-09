# -*- coding: utf-8 -*-

from StepBase import Step, Configure
import subprocess
import os

class MergeBamAlign(Step):
    def __init__(self,
                 unmappedBamInput = None,
                 alignedBamInput = None,
                 bamOutputDir = None,
                 refSequence = None,
                 secondAlign = False,
                 pairedRun = False,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('unmappedBamInput', unmappedBamInput)
        self.setParamIO('alignedBamInput', alignedBamInput)
        self.setParamIO('bamOutputDir', bamOutputDir)
        self.setParamIO('refSequence', refSequence)

        self.initIO()

        self.setParam('secondAlign', secondAlign)
        self.setParam('pairedRun', pairedRun)

    def impInitIO(self,):
        unmappedBamInput = self.getParamIO('unmappedBamInput')
        alignedBamInput = self.getParamIO('alignedBamInput')
        bamOutputDir = self.getParamIO('bamOutputDir')
        refSequence = self.getParamIO('refSequence')
        
        self.setInputDirOrFile('unmappedBamInput', unmappedBamInput)
        self.setInputDirOrFile('alignedBamInput', alignedBamInput)
        self.setInputDirOrFile('refSequence', refSequence)
        self.setOutputDirNTo1('bamOutput', os.path.join(bamOutputDir, 'merged.bam'), '', 'unmappedBamInput')

        if unmappedBamInput is not None:
            self._setInputSize(len(self.getInputList('unmappedBamInput')))

    def call(self, *args):
        unmappedBamUpstream = args[0]
        alignedBamUpstream = args[1]

        self.setParamIO('unmappedBamInput', unmappedBamUpstream.getOutput('bamOutput'))
        self.setParamIO('alignedBamInput', alignedBamUpstream.getOutput('bamOutput'))

    def _singleRun(self,i):
        unmappedBamInput = self.getInputList('unmappedBamInput')
        alignedBamInput = self.getInputList('alignedBamInput')
        bamOutput = self.getOutputList('bamOutput')
        refSequence = self.getInputList('refSequence')
        
        secondAlign = self.getParam('secondAlign')
        pairedRun = self.getParam('pairedRun')

        cmdline = [
                'java -Xmx4g -jar ../../dropseq/Drop-seq_tools-1.13/jar/lib/picard-2.10.3.jar MergeBamAlignment',
                'UNMAPPED_BAM=%s'%(unmappedBamInput[i]), 'ALIGNED_BAM=%s'%(alignedBamInput[i]), 'OUTPUT=%s'%(bamOutput[i]),
                'REFERENCE_SEQUENCE=%s'%(refSequence[i]), 'INCLUDE_SECONDARY_ALIGNMENTS=%s'%(str(secondAlign)),
                'PAIRED_RUN=%s'%(str(pairedRun))
        ]
        self.callCmdline(cmdline)
