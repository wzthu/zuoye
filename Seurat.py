# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from stepbase import Step,Configure
import subprocess
import os

class Seurat(Step):
    def __init__(self,
                 inputDir = None,                 
                 outputDir = None, 
                 threads = None,
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        """
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function.
        """
        
        # set all input and output parameters
        self.setParamIO('inputDir',inputDir)
        self.setParamIO('outputDir',outputDir)   
        # call self.initIO()
        self.initIO()
        #set other parameters
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)
        self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        inputDir = self.getParamIO('inputDir')
        outputDir = self.getParamIO('outputDir')        
        if inputDir is None:
        	self.setInput('barcodes',None)
        	self.setInput('genes', None)
        	self.setInput('matrix',None)
        else:
        	self.setInput('barcodes',os.path.join(inputDir,'barcodes.tsv'))
        	self.setInput('genes',os.path.join(inputDir,'genes.tsv'))
        	self.setInput('matrix',os.path.join(inputDir,'matrix.mtx'))
        #set all input files        
   
        # create output file paths and set
        self.setOutputDir1To1('VlnPlot', outputDir,'VlnPlot','jpg','genes') 
        self.setOutputDir1To1('GenePlot', outputDir,'GenePlot','jpg','genes') 
        '''
    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        fastqUpstream = args[0]      
        
        # set all required input parameters from upstream object
        self.setParamIO('fastqInput1',fastqUpstream.getOutput('fastqOutput1'))
        '''


    def _multiRun(self,):
        barcodes = self.getInput('barcodes')
        genes = self.getInput('genes')
        matrix = self.getInput('matrix')
        VlnPlot = self.getOutput('VlnPlot')
        GenePlot = self.getOutput('GenePlot')
        cmdline = ['Rscript',
        		   'Seurat.R',
        			self.getParamIO('inputDir'),
                    self.getParamIO('outputDir')
        			]
        # cmdline = ' '.join(cmdline)		     
        # print(' '.join(cmdline))      
        self.callCmdline(cmdline)
           
            
