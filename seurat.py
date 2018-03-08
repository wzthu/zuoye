#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from StepBase import Step
import os

class Seurat(Step):
    def __init__(self,
                 outputdir = None,
                 barcodes = None,
                 genes = None,
                 rscript = None,
                 matrix = None,
                 densematrix = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        # if densematrix or sparsematrix?
        if matrix!=None and barcodes!=None and genes!=None:
            self.setParamIO('matrix', matrix)
            self.setParamIO('barcodes', barcodes)
            self.setParamIO('genes', genes)
        else:
            self.setParamIO('densematrix', densematrix)
        self.setParamIO('outputdir', outputdir)
        self.setParamIO('rscript', rscript)
        self.initIO()
        self._setMultiRun()

    def impInitIO(self,):

        matrix = self.getParamIO('matrix')
        barcodes = self.setParamIO('barcodes')
        genes = self.setParamIO('genes')
        outputdir = self.getParamIO('outputdir')

        # set input paths
        self.setInputDirOrFile('matrix', matrix)
        self.setInputDirOrFile('barcodes', barcodes)
        self.setInputDirOrFile('genes', genes)

        # set output paths
        self.setOutputDirNTo1('violinplot', os.path.join(outputdir, 'violinplot.jpeg'), '', 'genes')
        self.setOutputDirNTo1('geneplot', os.path.join(outputdir, 'geneplot.jpeg'), '', 'genes')
        self.setOutputDirNTo1('Elbowplot', os.path.join(outputdir, 'Elbowplot.jpeg'), '', 'genes')

        # need to be segmented
        self.setOutputDirNTo1('TSNEplot', os.path.join(outputdir, 'TSNEplot.jpeg'), '', 'genes')


    def call(self, *args):

        # the first object
        fastqUpstream = args[0]

        # set all required input parameters from upstream object
        self.setParamIO('barcodes', fastqUpstream.getOutput('barcodes'))
        self.setParamIO('genes', fastqUpstream.getOutput('genes'))
        self.setParamIO('matrix', fastqUpstream.getOutput('matrix'))

        # print(self.cmdline)

    def _multiRun(self, ):
        # get input parameters
        barcodes = self.getParamIO('barcodes')
        genes = self.getParamIO('genes')
        matrix = self.getParamIO('matrix')
        rscript = self.getParamIO('rscript')

        # get output parameters
        outputdir = self.getParamIO('outputdir')

        cmdline = ['Rscript',
                   rscript,
                   matrix,
                   barcodes,
                   genes,
                   outputdir]

        print(cmdline)
        self.callCmdline(cmdline)

