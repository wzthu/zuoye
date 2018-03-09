#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from StepBase import Step
import os

class Cellranger(Step):
    def __init__(self,
                 outputdir = None,
                 fastqInput = None,
                 refile = None,
                 expectcells = 2000,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        self.setParamIO('fastqInput', fastqInput)
        self.setParamIO('refile', refile)
        self.setParamIO('outputdir', outputdir)
        self.setParam('id', outputdir)
        self.initIO()
		
		
        self.setParam('expectcells', expectcells)
        self._setMultiRun()

    def impInitIO(self,):        
        fastqInput = self.getParamIO('fastqInput')
        outputdir = self.getParamIO('outputdir')
        refile = self.getParamIO('refile')

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

        self.setOutputDirNTo1('barcodes', os.path.join(outputdir, 'outs', 'filtered_gene_bc_matrices', 'barcode.tsv'), '', 'fastqInput')
        self.setOutputDirNTo1('genes', os.path.join(outputdir, 'outs', 'filtered_gene_bc_matrices', 'genes.tsv'), '', 'fastqInput')
        self.setOutputDirNTo1('matrix', os.path.join(outputdir, 'outs', 'filtered_gene_bc_matrices', 'matrix.mtx'), '', 'fastqInput')

    def call(self, *args):
        pass
        # print(self.cmdline)

    def _multiRun(self, ):
        fastqInput = self.getParamIO('fastqInput')
        refile = self.getParamIO('refile')
        id = self.getParam('id')
        outputdir = self.getParamIO('outputdir')
        expectcells = self.getParam('expectcells')

        # if same outputdir already exist, delete it
        if os.path.isdir(outputdir):
            self.cmdline1 = ['rm', '-r %s' % outputdir]
            self.callCmdline(self.cmdline1)

        # run 10x cellranger
        self.cmdline2 = ['cellranger', 'count', '--id=%s --expect-cells=%s --transcriptome=%s --fastqs=%s'
                         % (id, str(expectcells), refile, fastqInput)]
        self.callCmdline(self.cmdline2)

