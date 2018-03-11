#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from StepBase import Step, Configure
import os

class Seurat(Step):
    def __init__(self,
                 outputdir=None,
                 # barcodes=None,
                 # genes=None,
                 inputdir = None,
                 rscript=None,
                 # matrix=None,
                 # densematrix = None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)
        # if densematrix or sparsematrix?
        # if matrix!=None and barcodes!=None and genes!=None:
        # self.setParamIO('matrix', matrix)
        # self.setParamIO('barcodes', barcodes)
        # self.setParamIO('genes', genes)
        # else:
        #     self.setParamIO('densematrix', densematrix)
        self.setParamIO('inputdir', inputdir)
        self.setParamIO('outputdir', outputdir)
        self.setParamIO('rscript', rscript)

        self.initIO()
        self._setMultiRun()

    def impInitIO(self, ):
        # matrix = self.getParamIO('matrix')
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        outputdir = self.getParamIO('outputdir')
        inputdir = self.getParamIO('inputdir')

        # if outputdir is None, os will error
        if outputdir is None:
           self.setParamIO('outputdir', Configure.getTmpDir())
           outputdir = self.getParamIO('outputdir')

        # set output/input paths
        # if inputdir is None, os will error
        if inputdir is not None:
            self.setInputDirOrFile('barcodes', os.path.join(inputdir, 'barcodes.tsv'))
            self.setInputDirOrFile('genes', os.path.join(inputdir, 'genes.tsv'))
            self.setInputDirOrFile('matrix', os.path.join(inputdir, 'matrix.mtx'))
            self.setOutputDirNTo1('violinplot', os.path.join(outputdir, 'violinplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('geneplot', os.path.join(outputdir, 'geneplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('variableGenes', os.path.join(outputdir, 'variableGenes.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('Elbowplot', os.path.join(outputdir, 'Elbowplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('TSNEplot', os.path.join(outputdir, 'TSNEplot.jpeg'), '', 'barcodes')

    def call(self, *args):
        # the first object
        cellrangerUpstream = args[0]
        # set all required input parameters from upstream objects
        self.setParamIO('inputdir', cellrangerUpstream.getOutput('finaldir'))
        # self.setParamIO('genes', cellrangerUpstream.getOutput('genes'))
        # self.setParamIO('matrix', cellrangerUpstream.getOutput('matrix'))

        # print(self.cmdline)

    def _multiRun(self, ):
        # get input parameters
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        # matrix = self.getParamIO('matrix')
        inputdir = self.getParamIO('inputdir')
        rscript = self.getParam('rscript')
        # get output parameters
        outputdir = self.getParamIO('outputdir')

        # replace file path
        inputdir = inputdir.replace('/home/cfeng', '/data')
        # matrix = matrix.replace('/home/cfeng', '/data')
        rscript = rscript.replace('/home/cfeng', '/data')
        # barcodes = barcodes.replace('/home/cfeng', '/data')

        outputdir = outputdir.replace('/home/cfeng', '/data')

        cmdline = ['docker',
                   'exec',
                   '-it',
                   'hca_fengchen',
                   '/root/software/anaconda3/bin/Rscript',
                   rscript,
                   inputdir,
                   outputdir]
        print(cmdline)
        self.callCmdline(cmdline)

