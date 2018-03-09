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
                 # densematrix = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)
        # if densematrix or sparsematrix?
        # if matrix!=None and barcodes!=None and genes!=None:
        self.setParamIO('matrix', matrix)
        self.setParamIO('barcodes', barcodes)
        self.setParamIO('genes', genes)
        # else:
        #     self.setParamIO('densematrix', densematrix)
        self.setParamIO('outputdir', outputdir)
        self.setParamIO('rscript', rscript)
        self.initIO()
        self._setMultiRun()

    def impInitIO(self,):

        matrix = self.getParamIO('matrix')
        barcodes = self.getParamIO('barcodes')
        genes = self.getParamIO('genes')
        outputdir = self.getParamIO('outputdir')

        # set input paths
        self.setInputDirOrFile('matrix', matrix)
        self.setInputDirOrFile('barcodes', barcodes)
        self.setInputDirOrFile('genes', genes)

        # set output paths
        self.setOutputDir1To1('violinplot', outputdir, 'violinplot', 'jpeg', 'genes')
        self.setOutputDir1To1('geneplot', outputdir, 'geneplot', 'jpeg', 'genes')
        self.setOutputDir1To1('Elbowplot', outputdir, 'Elbowplot', 'jpeg', 'genes')
        # need to be segmented
        self.setOutputDir1To1('TSNEplot', outputdir, 'TSNEplot', 'jpeg', 'genes')


    def call(self, *args):

        # the first object
        cellrangerUpstream = args[0]

        # set all required input parameters from upstream object
        self.setParamIO('barcodes', cellrangerUpstream.getOutput('barcodes'))
        self.setParamIO('genes', cellrangerUpstream.getOutput('genes'))
        self.setParamIO('matrix', cellrangerUpstream.getOutput('matrix'))

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

