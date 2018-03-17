#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from StepBase import Step, Configure
from Cellranger import Cellranger
import os

class Seuratrun(Step):
    def __init__(self,
                 outputdir=None,
                 # barcodes=None,
                 # genes=None,
                 inputdir = None,
                 rscript=None,
                 # datatype = None,
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
        # self.setParam('datatype', datatype)

        self.initIO()
        self._setMultiRun()

    def impInitIO(self, ):
        # matrix = self.getParamIO('matrix')
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        outputdir = self.getParamIO('outputdir')
        inputdir = self.getParamIO('inputdir')
        rscript = self.getParamIO('rscript')

        self.setInputDirOrFile('rscript', rscript)
        # if outputdir is None, os will error
        if outputdir is None:
           self.setParamIO('outputdir', Configure.getTmpDir())
           outputdir = self.getParamIO('outputdir')

        # set output/input paths
        # if inputdir is None, os will error
        if inputdir is not None:
            self.setInputDirOrFile('inputfile', os.path.join(inputdir, 'Preprocessing.Rdata'))
            # self.setInputDirOrFile('barcodes', os.path.join(inputdir, 'barcodes.tsv'))
            # self.setInputDirOrFile('genes', os.path.join(inputdir, 'genes.tsv'))
            # self.setInputDirOrFile('matrix', os.path.join(inputdir, 'matrix.mtx'))
            # self.setOutputDirNTo1('violinplot', os.path.join(outputdir, 'violinplot.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('geneplot', os.path.join(outputdir, 'geneplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('variableGenes', os.path.join(outputdir, 'variableGenes.jpeg'), '', 'inputfile')
            self.setOutputDirNTo1('Elbowplot', os.path.join(outputdir, 'Elbowplot.jpeg'), '', 'inputfile')
            self.setOutputDirNTo1('TSNEplot', os.path.join(outputdir, 'TSNEplot.jpeg'), '', 'inputfile')

    def call(self, *args):
        # the first object
        cellrangerUpstream = args[0]
        # set all required input parameters from upstream objects
        # if isinstance(cellrangerUpstream, Cellranger):
        #     datatype = '10x'
        #     self.setParam('datatype', datatype)
        self.setParamIO('inputdir', cellrangerUpstream.getParamIO('outputdir'))
        # self.setParamIO('genes', cellrangerUpstream.getOutput('genes'))
        # self.setParamIO('matrix', cellrangerUpstream.getOutput('matrix'))

        # print(self.cmdline)

    def _multiRun(self, ):
        # get input parameters
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        # matrix = self.getParamIO('matrix')
        inputfile = self.getInput('inputfile')
        rscript = self.getParamIO('rscript')
        # datatype = self.getParam('datatype')
        # get output parameters
        outputdir = self.getParamIO('outputdir')

        cmdline = ['Rscript',
                   rscript,
                   inputfile,
                   # datatype,
                   outputdir]
        print(cmdline)
        self.callCmdline(cmdline=cmdline, dockerVersion='V1')
        #
        # replace file path
        # inputdir = inputdir.replace('/home/cfeng', '/data')
        # # matrix = matrix.replace('/home/cfeng', '/data')
        # rscript = rscript.replace('/home/cfeng', '/data')
        # # barcodes = barcodes.replace('/home/cfeng', '/data')
        #
        # outputdir = outputdir.replace('/home/cfeng', '/data')
        # if datatype is ('10x' or '10X' or 'tenx' or 'Tenx' or 'TenX' or 'tenX'):
        #     cmdline = ['Rscript',
        #                rscript,
        #                inputdir,
        #                #datatype,
        #                outputdir]
        #     print(cmdline)
        #     self.callCmdline(cmdline=cmdline, dockerVersion='V1')
        # else:
        #     cmdline = ['Rscript',
        #                rscript,
        #                inputdir,
        #                # datatype,
        #                outputdir]
        #     print(cmdline)
        #     self.callCmdline(cmdline=cmdline, dockerVersion='V1')

    def getMarkdownEN(self):
        rmd = '''
---
title: "Seuratrun"
author: "Author:Frankie"
output: html_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
`Runobject = Seuratrun(inputdir, outputdir, rscript)(upstream)`

### Input: 
1. Created from original data
- inputdir : Path of upstream files(Seurat object).
- outputdir: Path of all the analysis results
- rscript : Path of alternative Rscript.
2. From upstream object such as Seuratpreprocessing object
- rscript : Path of alternative Rscript.

### Output:
- All the results will be saved in the outputdir
- tSne plot& Variablegenes plot& PCAelbowplot

### Variablegenes
![]({variableimagepath})

### PCAelbowplot
![]({elbowimagepath})

### tSne plot
![]({tsneimagepath})
        '''.format(variableimagepath=self.getOutput('variableGenes'),
                   elbowimagepath=self.getOutput('Elbowplot'),
                   tsneimagepath=self.getOutput('TSNEplot'))
        return rmd
