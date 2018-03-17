#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from StepBase import Step, Configure
from Cellranger import Cellranger
import os

class Seuratpreprocessing(Step):
    def __init__(self,
                 outputdir=None,
                 # barcodes=None,
                 # genes=None,
                 inputdir = None,
                 rscript=None,
                 datatype = None,
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
        self.setParam('datatype', datatype)

        self.initIO()
        self._setMultiRun()

    def impInitIO(self, ):
        # matrix = self.getParamIO('matrix')
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        outputdir = self.getParamIO('outputdir')
        inputdir = self.getParamIO('inputdir')
        rscript = self.getParamIO('rscript')
        datatype = self.getParam('datatype')

        self.setInputDirOrFile('rscript', rscript)
        # if outputdir is None, os will error
        if outputdir is None:
           self.setParamIO('outputdir', Configure.getTmpDir())
           outputdir = self.getParamIO('outputdir')

        # set output/input paths
        # if inputdir is None, os will error
        if inputdir is not None:
            if datatype is'10x':
                self.setInputDirOrFile('barcodes', os.path.join(inputdir, 'barcodes.tsv'))
                self.setInputDirOrFile('genes', os.path.join(inputdir, 'genes.tsv'))
                self.setInputDirOrFile('matrix', os.path.join(inputdir, 'matrix.mtx'))
            else:
                self.setInputDirOrFile('matrix', inputdir)

            self.setOutputDirNTo1('Rdata', os.path.join(outputdir, 'Preprocessing.Rdata'), '', 'barcodes')
            self.setOutputDirNTo1('violinplot', os.path.join(outputdir, 'violinplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('geneplot', os.path.join(outputdir, 'geneplot.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('variableGenes', os.path.join(outputdir, 'variableGenes.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('Elbowplot', os.path.join(outputdir, 'Elbowplot.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('TSNEplot', os.path.join(outputdir, 'TSNEplot.jpeg'), '', 'barcodes')

    def call(self, *args):
        # the first object
        Upstream = args[0]
        # set all required input parameters from upstream objects
        if isinstance(Upstream, Cellranger):
            datatype = '10x'
            self.setParam('datatype', datatype)
        else:
            self.setParam('datatype', 'densemtx')
        self.setParamIO('inputdir', Upstream.getParamIO('finaldir'))
        # self.setParamIO('genes', cellrangerUpstream.getOutput('genes'))
        # self.setParamIO('matrix', cellrangerUpstream.getOutput('matrix'))

        # print(self.cmdline)

    def _multiRun(self, ):
        # get input parameters
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        # matrix = self.getParamIO('matrix')
        inputdir = self.getParamIO('inputdir')
        rscript = self.getParamIO('rscript')
        datatype = self.getParam('datatype')
        # get output parameters
        outputdir = self.getParamIO('outputdir')

        #
        # replace file path
        # inputdir = inputdir.replace('/home/cfeng', '/data')
        # # matrix = matrix.replace('/home/cfeng', '/data')
        # rscript = rscript.replace('/home/cfeng', '/data')
        # # barcodes = barcodes.replace('/home/cfeng', '/data')
        #
        # outputdir = outputdir.replace('/home/cfeng', '/data')
        cmdline = ['Rscript',
                   rscript,
                   inputdir,
                   datatype,
                   outputdir]
        print(cmdline)
        self.callCmdline(cmdline=cmdline, dockerVersion='V1')

    def getMarkdownEN(self):
        rmd = '''
---
title: "Seuratpreprocessing"
author: "Author:Frankie"
output: html_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
`Preobject = Seuratpreprocessing(inputdir, outputdir, datatype=['10x','densemtx'], rscript)(upstream)`

### Input: 
1. Created from original data
- inputdir : Path of Matrix directory(one dense matrix or sparse matrix per time).
- outputdir: Path of all the analysis results
- datatpe: Specify the inputdata type, such as '10x'
- rscript : Path of alternative Rscript.
2. From upstream object such as Seuratpreprocessing object
- rscript : Path of alternative Rscript.

### Output:
- Seurat object saved in outputdir for downstream analysis
- Violinplot and geneplot saved in outputdir for you to do filtering

### Violinplot
![]({vioimagepath})

### Geneplot
![]({geneimagepath})
        '''.format(vioimagepath=self.getOutput('violinplot'), geneimagepath = self.getOutput('geneplot'))
        return rmd

