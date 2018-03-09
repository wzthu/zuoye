# -*- coding: utf-8 -*-
"""
Created on Tus Mar  6 12:28:52 2018

@author: ShengquanChen
"""

from StepBase import Step,Configure
import subprocess
import os

class FastqDump(Step):
    def __init__(self,
                 sraInput1 = None,              
                 fastqOutputDir = None, 
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        
        self.setParamIO('sraInput1',sraInput1)
        self.setParamIO('fastqOutputDir',fastqOutputDir)

       
        self.initIO()
            
        
    def impInitIO(self,):        
        sraInput1 = self.getParamIO('sraInput1')
        fastqOutputDir = self.getParamIO('fastqOutputDir')

        #set all input files        
        self.setInputDirOrFile('sraInput1',sraInput1) 
       
        # self.setOutputDir1To1('fastqOutputDir', fastqOutputDir,'fastqDump','fastq','sraInput1',sep='_') 
        self.setOutputDir1To1('fastqOutput1', fastqOutputDir, None, '1.fastq','sraInput1',sep='_') 
        self.setOutputDir1To1('fastqOutput2', fastqOutputDir, None, '2.fastq','sraInput1',sep='_') 
        
        if sraInput1 is not None:
            self._setInputSize(len(self.getInputList('sraInput1')))
        
    def call(self, *args):
        fastqUpstream = args[0]      
        
        self.setParamIO('sraInput1',fastqUpstream.getOutput('sraOutput1'))
            
    def _singleRun(self, i):
        sraInput1 = self.getInputList('sraInput1')
        fastqOutputDir = self.getParamIO('fastqOutputDir')
        cmdline = ['fastq-dump', '--split-3',
                    sraInput1[i],
                    '-O', fastqOutputDir
                    ]
                    
        result = self.callCmdline(cmdline)
        # f = open(mapRsOutput[i],'wb')   
        # f.write(result.stderr)
            
        