#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
20180306
"""
from stepbase import Step,Configure
import subprocess
import os

class Cellranger(Step):
    def __init__(self,
        outputdir = None,
        fastqInput = None, 
        refile = None,
        expectcells = 2000,
        cmdParam = None,
        #threads = None,
        **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        self.setParamIO('fastqInput',fastqInput)
        self.setParamIO('refile',refile)
        self.setParamIO('outputdir',outputdir)
        self.initIO()
        # self.setPaif threads is None:
        # threads = Configure.getThreads()
        # self.setParam('threads',threads)
        self.setParam('expectcells', expectcells)
        self._setMultiRun()

	#  `

    def impInitIO(self,):        
        fastqInput = self.getParamIO('fastqInput')
        #fastqInput2 = self.getParamIO('fastqInput2')
        refile = self.getParamIO('refile')
        outputdir = self.getParamIO('outputdir')
        
        #set all input files
        
        self.setInputDirOrFile('fastqInput',fastqInput)
        #self.setInputDir('fastqInput',fastqInput)
        self.setInputDirOrFile('refile',refile)  
        
        self.setOutDirNTo1('barcodes',os.path.join(outputdir,'outs','filtered_gene_bc_matrices',Configure.getGenome(),'barcode.tsv','fastqInput')
        self.setOutDirNTo1('genes',os.path.join(outputdir,'outs','filtered_gene_bc_matrices',Configure.getGenome(),'genes.tsv','fastqInput')
        self.setOutDirNTo1('matrix',os.path.join(outputdir,'outs','filtered_gene_bc_matrices',Configure.getGenome(),'matrix.mtx','fastqInput')
        
            
#    def setOutputDir1To1(self, outputName, outputDir, outputPrefix, outputSuffix, inputName):
#        # inputList = self.getInput(inputName)
#        # if inputList is None:
#            # self.setOutput(outputName, None)
#        # else:
#            # if not isinstance(inputList,list):
#                # inputList = [inputList]
#            # outputList = []            
#            # for i in range(len(inputList)):
#                # outputList.append(outputPrefix + '_' + str(i) + '_' + outputSuffix)
#            # if outputDir is None:
#                # self.setOutput(outputName,Configure.getTmpPath(outputList))
#            # else:        
#				self.setOutput(outputName,[outputPrefix + '_' + str(i) + '_' + outputSuffix])
				
#    def _multiRun(self,):
#        fastqInput = self.getInputList('fastqInput')
#        expectcells = self.getParam('expectcells')
#        refile = self.getInputList('refile')
#        outputdir = self.getOutputList('outputdir')
#        cmdline = ['cellranger count',
#           '--id=%s --expect-cells=%s --transcriptome=%s --fastqs=%s' 
#           % (outputdir,str(expectcells),refile,fastqInput)]
#        print(cmdline)
#        self.callCmdline(cmdline)
        
    def call(self, *args):
        pass #print(self.cmdline)          
	
    def _singleRun(self, i):
        fastqInput = self.getInput('fastqInput')
        expectcells = self.getParam('expectcells')
        refile = self.getInput('refile')#getParamIO('refile')
        outpout = self.getParamIO('refile')
        self.cmdline = ["cellranger count --id=%s --expect-cells=%s --transcriptome=%s --fastqs=%s" % (output[i],str(expectcells),refile[i],fastqInput[i])]
        print(self.cmdline[0])
        self.callCmdline(self.cmdline[0])