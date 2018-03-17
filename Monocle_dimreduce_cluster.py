# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from StepBase import Step,Configure
import subprocess
import os

class Monocle_dimreduce_cluster(Step):
    def __init__(self,
                 imageRdata = None,
                 outputpath = None,
                 num_PCA = 3,
                 cluster_num = 2,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        """
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function.
        """
        
        #Configure.enableDocker(False)
        # set all input and output parameters
        self.setParamIO('imageRdata',imageRdata)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters

        self.setParam('num_PCA',num_PCA)
        self.setParam('cluster_num',cluster_num)
        self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        imageRdata = self.getParamIO('imageRdata')        
        outputpath = self.getParamIO('outputpath')  

        
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  
        self.setInput('imageRdata', imageRdata)
        # create output file paths and set
        self.setOutputDir1To1('densitypeak_cluster',outputpath, 'densitypeak_cluster','jpg','imageRdata')
       
    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        MonocleQCupstream = args[0]      
        
        # set all required input parameters from upstream object
        self.setParamIO('imageRdata',MonocleQCupstream.getOutput('MonocleQCimage'))
        #print(MonocleQCupstream.getOutput('MonocleQCimage'))
        
    def _multiRun(self,):
        imageRdata = self.getInput('imageRdata')
        num_PCA = self.getParam('num_PCA')
        cluster_num = self.getParam('cluster_num')
        cmdline = ['Rscript',
                    '/data/Monocle_dimreduce_cluster.R',
        			imageRdata[0],
                    str(num_PCA),
                    str(cluster_num),
                    self.getParamIO('outputpath')
        			]
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)