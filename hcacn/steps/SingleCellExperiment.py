# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from ..core import Step,Configure

import os

class SingleCellExperiment(Step):
    def __init__(self,
                 matrix_file = None,
                 ann_file = None,
                 matrix_format = None,
                 outputpath = None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        
         

        # set all input and output parameters
        self.setInput('matrix_file',matrix_file)
        self.setInput('ann_file',ann_file)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters
        self.setParam("matrix_format",matrix_format)
    def impInitIO(self,):
        # obtain all input and output parameters           
        outputpath = self.getParamIO('outputpath')  
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath') 

        matrix_file = self.getInput('matrix_file') 
        # create output file paths and set
        self.setOutputDir1To1('sceOutput',outputpath,None,".RData","matrix_file",sep='')

        # Rscripts
        self.setInputRscript('Rscript','SingleCellExperiment.R')
        if matrix_file is not None:
            self._setInputSize(len(self.getInputList('matrix_file')))

    def call(self,*args):

        Upstream = args[0]
        '''
        No Implement till now!
        '''
        #other things

    def getMarkdownEN(self,):
        mdtext = """
## SC3_DE Usage

SC3_DE('/path/to/matrix.csv','/path/to/annotation.csv','matrix_format','/path/to/outputDir')  

## SingleCellExperiment Result  
Successful generate SingleCellExperiment from matrix.
"""

        return mdtext
            
            
    def _singleRun(self,i):
        # obtain all input and output dir list
        matrix_files = self.getInputList('matrix_file')
        ann_files = self.getInputList('ann_file')
        sceOutputDirs = self.getOutputDir('sceOutput')
        matrix_format = self.getParam('matrix_format')
        assert matrix_format in ['LOG','ORIGIN']
        Rscript = self.getInput('Rscript')
        cmdline =['Rscript',
                  Rscript,
                   matrix_files[i],
                   ann_files[i],
                  matrix_format,
                   sceOutputDirs[i]
                   ]
        self.callCmdline('V1', cmdline)