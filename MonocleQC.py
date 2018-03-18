# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from StepBase import Step,Configure
import subprocess
import os

class MonocleQC(Step):
    def __init__(self,
                 matrixdata = None,
                 outputpath = None,
                 min_expression = 0.1,
                 num_cells_expressed_threshold = 5,
                 TotalmRNAs = 1e6, 
                 mean_expression_threshold=0.1,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        """
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function.
        """
        
        #Configure.enableDocker(False)
        # set all input and output parameters
        self.setParamIO('matrixdata',matrixdata)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters
        #if threads is None:
        #    threads = Configure.getThreads()
        self.setParam('min_expression',min_expression)
        self.setParam('num_cells_expressed_threshold',num_cells_expressed_threshold)
        self.setParam('TotalmRNAs',TotalmRNAs)
        self.setParam('mean_expression_threshold',mean_expression_threshold)
        #if threads is None:
        #    threads=Configure.getThreads()
        #self.setParam('threads', threads)
        self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        matrixdata = self.getParamIO('matrixdata')        
        outputpath = self.getParamIO('outputpath')  

        
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  
        self.setInput('matrixdata', matrixdata)
        # create output file paths and set
        #if outputpath is None:
        #   self.setOutputDir1To1('density_Total_mRNAs',None)
        #else:
        #   self.setOutputDir1To1('samOutput', samOutputDir,None,'sam','fastqInput1') 
        #self.setOutputDirNTo1('density_Total_mRNAs', None, Configure.getTmpPath('density_Total_mRNAs.jpg'),'matrixdata')
        #self.setOutputDirNTo1('meanexpression_disersionemprical', None, Configure.getTmpPath('meanexpression_disersionemprical.jpg'),'matrixdata')
        # self.setOutputDirNTo1('PCvariance', None, Configure.getTmpPath('PCvariance.jpg'),'matrixdata')
        self.setOutputDir1To1('density_Total_mRNAs',outputpath, 'density_Total_mRNAs','jpg','matrixdata')
        self.setOutputDir1To1('meanexpression_disersionemprical',outputpath,'meanexpression_disersionemprical','jpg','matrixdata')
        self.setOutputDir1To1('PCvariance', outputpath,'PCvariance','jpg','matrixdata')
        self.setOutputDir1To1('MonocleQCimage', outputpath,'MonocleQCimage','Rdata','matrixdata')
    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        dropseqUpstream = args[0]      
        
        # set all required input parameters from upstream object
        self.setParamIO('matrixdata',dropseqUpstream.getOutput('dgeOutput'))
        
    def _multiRun(self,):
        matrixdata = self.getInput('matrixdata')
        min_expression = self.getParam('min_expression')
        num_cells_expressed_threshold = self.getParam('num_cells_expressed_threshold')
        TotalmRNAs = self.getParam('TotalmRNAs')
        mean_expression_threshold = self.getParam('mean_expression_threshold')
        cmdline = ['Rscript',
                    '/data/MonocleQC.R',
        			matrixdata,
                    str(min_expression),
                    str(num_cells_expressed_threshold),
                    str(TotalmRNAs),
                    str(mean_expression_threshold),
                    self.getParamIO('outputpath')
        			]
        '''cmdline='Rscript 
        /dataMonocleQC.R %s %d %d %d %d %s'%(matrixdata,
        min_expression,num_cells_expressed_threshold,
        TotalmRNAs,
        mean_expression_threshold,
        self.getParamIO('outputpath'))'''
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)

    def getMarkdownEn(self,):
        mdtext="""
        ### Monocle QC Result
        The Total_mRNAs~density Curve of all cells is shown below:
        ```{{r setup, include=FALSE}}
        knitr::opts_chunk$set(echo = TRUE)
        ```
        ### Monocle QC Result
        ```{{r,eval=FALSE}}
        #don't run
        MonocleQC(matrixdata='', outputpath='') 
        ```
        ####The Total_mRNAs~density Curve of all cells is shown below:

        ![]({density_Total_mRNAs})

        ####The picture below shows how variability (dispersion) in a gene's expression depends on the average expression across cells:
        ![]({meanexpression_disersionemprical})

        ####The picture below shows the variance explained by each component:
        ![]({PCvariance})

        """.format(density_Total_mRNAs = self.getOutput('density_Total_mRNAs'),
                   meanexpression_disersionemprical = self.getOutput('meanexpression_disersionemprical'),
                   PCvariance = self.getOutput('PCvariance'))