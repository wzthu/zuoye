#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from StepBase import Step, Configure
import os

class Cellranger(Step):
    def __init__(self,
                 outputdir = None,
                 fastqInput = None,
                 refile = None,
                 expectcells = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        self.setParamIO('fastqInput', fastqInput)
        self.setParamIO('refile', refile)
        self.setParamIO('outputdir', outputdir)
        self.initIO()

        self.setParam('expectcells', expectcells)
        #Configure.enableDocker(False)
        self._setMultiRun()

    def impInitIO(self,):        
        fastqInput = self.getParamIO('fastqInput')
        refile = self.getParamIO('refile')
        outputdir = self.getParamIO('outputdir')

        self.setInputDirOrFile('fastqInput', fastqInput)
        self.setInputDirOrFile('version', os.path.join(refile, 'version'))
        self.setInputDirOrFile('Reference', os.path.join(refile, 'reference.json'))
        self.setInputDirOrFile('README', os.path.join(refile, 'README.BEFORE.MODIFYING'))
        for i in ['chrLength.txt', 'chrName.txt',
                  'exonGeTrInfo.tab', 'geneInfo.tab',
                  'genomeParameters.txt', 'SAindex',
                  'sjdbList.fromGTF.out.tab', 'transcriptInfo.tab',
                  'chrNameLength.txt', 'chrStart.txt',
                  'exonInfo.tab', 'Genome', 'SA',
                  'sjdbInfo.txt', 'sjdbList.out.tab']:
            self.setInputDirOrFile(i, os.path.join(refile, 'star', i))
        self.setInputDirOrFile('genes.pickle', os.path.join(refile, 'pickle', 'genes.pickle'))
        self.setInputDirOrFile('genes.gtf', os.path.join(refile, 'genes', 'genes.gtf'))
        self.setInputDirOrFile('genome.fa', os.path.join(refile, 'fasta', 'genome.fa'))

        if outputdir is None:
            self.setParamIO('outputdir', Configure.getTmpDir())
            outputdir = self.getParamIO('outputdir')
            self.resultdir = 'Cellranger'
        else:
            self.resultdir = ''

        self.setParamIO('finaldir', os.path.join(outputdir, self.resultdir, 'outs', 'filtered_gene_bc_matrices', 'hg19'))
        self.setOutputDirNTo1('genes', os.path.join(outputdir, self.resultdir,'outs', 'filtered_gene_bc_matrices', 'hg19', 'genes.tsv'), '', 'fastqInput')
        self.setOutputDirNTo1('matrix', os.path.join(outputdir, self.resultdir,'outs', 'filtered_gene_bc_matrices', 'hg19', 'matrix.mtx'), '', 'fastqInput')
        self.setOutputDirNTo1('barcodes', os.path.join(outputdir,self.resultdir, 'outs', 'filtered_gene_bc_matrices', 'hg19', 'barcodes.tsv'), '', 'fastqInput')

    def call(self, *args):
        pass
        # print(self.cmdline)

    def _multiRun(self, ):
        fastqInput = self.getParamIO('fastqInput')
        refile = self.getParamIO('refile')
        # id = self.getParam('id')
        outputdir = self.getParamIO('outputdir')
        expectcells = self.getParam('expectcells')
        # cd outputdir to run the command
        # self.cmdline1 = ['cd', '%s' % outputdir]
        # self.callCmdline(self.cmdline1)
        if expectcells is not None:
            if self.resultdir is '':
                # use given path
                dir = outputdir.split('/')[-1]
                if os.path.isdir(outputdir):
                    cmdline = ['rm -r', '%s' % outputdir]
                    self.callCmdline(cmdline=cmdline, dockerVersion='V1',stdoutToLog=False)
                cmdline = ['cellranger', 'count','--id=%s --expect-cells=%s --transcriptome=%s --fastqs=%s'
                           % (dir, str(expectcells), refile, fastqInput)]
                print(cmdline)
                self.callCmdline(cmdline= cmdline, dockerVersion='V1', stdoutToLog=False)
            else:
                # use default filepath
                cmdline = ['cd', outputdir, '&&', 'rm -rf Cellranger', '&&', 'cellranger', 'count',
                           '--id=Cellranger', '--expect-cells=%s --transcriptome=%s --fastqs=%s'
                           % (str(expectcells), refile, fastqInput)]
                self.callCmdline(cmdline=cmdline, dockerVersion='V1', stdoutToLog=False)

        else:
            if self.resultdir is '':
                # use given path
                dir = outputdir.split('/')[-1]
                if os.path.isdir(outputdir):
                    cmdline = ['rm -r', '%s' % outputdir]
                    self.callCmdline(cmdline=cmdline, stdoutToLog=False)
                cmdline = ['cellranger', 'count', '--id=%s  --transcriptome=%s --fastqs=%s'
                           % (dir, refile, fastqInput)]
                print(cmdline)
                self.callCmdline(cmdline=cmdline, dockerVersion='V1', stdoutToLog=False)
            else:
                # use default filepath
                # cd target dir
                # run 10x cellranger
                cmdline = ['cd', outputdir,'&&','rm -rf Cellranger', '&&', 'cellranger', 'count',
                           '--id=Cellranger', '--transcriptome=%s --fastqs=%s'
                           % (refile, fastqInput)]
                self.callCmdline(cmdline=cmdline, dockerVersion='V1', stdoutToLog=False)

    def getMarkdownEN(self):
        rmd = '''
---
title: "Cellranger"
author: "Author:Frankie"
output: html_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
`Cellrangerobject = Cellranger(fastqinput, outputdir, refile, expectcells)`

### Input: 
- fastqInput : Path of folder containing fastq files.
- outputdir: Path of all the analysis results
- refile : Path of folder containing 10x-compatible reference.
- [expect-cells] : Extra options that can tune the program(often as default)

### Output:
- All the results will be saved in the outputdir
- Expression matrix can be found in outputdir/outs/filtered_gene_bc_matrices
- Summary html can be found in 
        '''
        return rmd

