# -*- coding: utf-8 -*-
"""
Created on Tus Mar  7 12:28:52 2018

@author: ShengquanChen
"""

from StepBase import Step,Configure
import subprocess
import os

class HTSeq_sam2count(Step):
    def __init__(self,
                 samInput1 = None,  
                 gtfInput1 = None,             
                 countOutputDir = None, 
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        
        self.setParamIO('samInput1',samInput1)
        self.setParamIO('gtfInput1',gtfInput1)
        self.setParamIO('countOutputDir',countOutputDir)

       
        self.initIO()
            
        
    def impInitIO(self,):        
        samInput1 = self.getParamIO('samInput1')
        gtfInput1 = self.getParamIO('gtfInput1')
        countOutputDir = self.getParamIO('countOutputDir')
        if countOutputDir is None:
            self.setParamIO('countOutputDir',Configure.getTmpDir())
            

        #set all input files        
        self.setInputDirOrFile('samInput1',samInput1) 
        self.setInputDirOrFile('gtfInput1',gtfInput1) 
       
        self.setOutputDir1To1('countOutput', countOutputDir, None, 'count','samInput1') 
        
        if samInput1 is not None:
            self._setInputSize(len(self.getInputList('samInput1')))
        
    def call(self, *args):
        htseqUpstream = args[0]      
        
        self.setParamIO('samInput1',htseqUpstream.getOutput('samOutput1'))
        self.setParamIO('gtfInput1',htseqUpstream.getOutput('gtfOutput1'))
            
    def _singleRun(self, i):
        samInput1 = self.getInputList('samInput1')
        gtfInput1 = self.getInputList('gtfInput1')
        countOutput = self.getOutputList('countOutput')
        cmdline = ['python -m HTSeq.scripts.count', '-s no',
                    samInput1[i],
                    gtfInput1[i],
                    '>',
                    countOutput[i]
                    ]
                    
        result = self.callCmdline('V1', cmdline)
        # f = open(mapRsOutput[i],'wb')   
        # f.write(result.stderr)
            
        