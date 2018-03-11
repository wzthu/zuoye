# -*- coding: utf-8 -*-
"""
Created on Tus Mar  10 10:28:52 2018

@author: ShengquanChen
"""

from StepBase import Step,Configure
import subprocess
import os

class Cuffmerge(Step):
    def __init__(self,
                 faInput1 = None,  
                 gtfInput1 = None,  
                 assembliesInput1 = None,
                 threads = None,
                 gtfOutputDir = None, 
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        
        self.setParamIO('faInput1',faInput1)
        self.setParamIO('gtfInput1',gtfInput1)
        self.setParamIO('assembliesInput1',assembliesInput1)
        self.setParamIO('gtfOutputDir',gtfOutputDir)

       
        self.initIO()
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)
            
        
    def impInitIO(self,):        
        faInput1 = self.getParamIO('faInput1')
        gtfInput1 = self.getParamIO('gtfInput1')
        assembliesInput1 = self.getParamIO('assembliesInput1')
        gtfOutputDir = self.getParamIO('gtfOutputDir')

        #set all input files        
        self.setInputDirOrFile('assembliesInput1',assembliesInput1) 
       
        self.setOutputDir1To1('gtfOutputDir', gtfOutputDir, None, 'gtf','assembliesInput1') 
        
        if assembliesInput1 is not None:
            self._setInputSize(len(self.getInputList('assembliesInput1')))
        
    def call(self, *args):
        htseqUpstream = args[0]              
        self.setParamIO('assembliesInput1',htseqUpstream.getOutput('assembliesOutput1'))
        # self.setParamIO('gtfInput1',htseqUpstream.getOutput('gtfOutput1'))
            
    def _singleRun(self, i):
        faInput1 = self.getParamIO('faInput1')
        gtfInput1 = self.getParamIO('gtfInput1')
        assembliesInput1 = self.getInputList('assembliesInput1')
        gtfOutputDir = self.getOutputList('gtfOutputDir')
        cmdline = ['docker run --rm -v /home/hca/Docker/Common_data:/data hca:py2 cuffmerge',
                    '-g', gtfInput1,
                    '-s', faInput1,
                    '-o', gtfOutputDir[i],
                    '-p', str(self.getParam('threads')),
                    assembliesInput1[i]
                    ]
                    
        result = self.callCmdline(cmdline)
        # f = open(mapRsOutput[i],'wb')   
        # f.write(result.stderr)
        # 
        # docker run --rm -v /home/hca/Docker/Common_data:/data hca:py2 cuffmerge -g /data/sqchen/genome.gtf -s /data/sqchen/hg19.fa -o /data/sqchen -p 8 /data/sqchen/assemblies.txt
            
        