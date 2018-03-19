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
                 assembliesInput = None,
                 threads = None,
                 gtfOutputDir = None, 
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        
        self.setParamIO('faInput1',faInput1)
        self.setParamIO('gtfInput1',gtfInput1)
        self.setParamIO('assembliesInput',assembliesInput)
        self.setParamIO('gtfOutputDir',gtfOutputDir)

       
        self.initIO()
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)
            
        
    def impInitIO(self,):        
        faInput1 = self.getParamIO('faInput1')
        gtfInput1 = self.getParamIO('gtfInput1')
        assembliesInput = self.getParamIO('assembliesInput')
        gtfOutputDir = self.getParamIO('gtfOutputDir')
        if gtfOutputDir is None:
            self.setParamIO('gtfOutputDir',Configure.getTmpDir())

        #set all input files        
        self.setInputDirOrFile('assembliesInput',assembliesInput) 
        
        if faInput1 is None:
            faInput1=Configure.getConfig('')
            self.setIput('faInput1',faInput1)
            self.setParamIO('faInput1',faInput1)
        else:
            self.setInput('faInput1',faInput1)

        if gtfInput1 is None:
            gtfInput1=Configure.getConfig('')
            self.setIput('gtfInput1',gtfInput1)
            self.setParamIO('gtfInput1',gtfInput1)
        else:
            self.setInput('gtfInput1',gtfInput1)

        if assembliesInput is not None:
            self._setInputSize(len(self.getInputList('assembliesInput')))
            merged_gtf=list()
            for i in range(len(self.getInputList('assembliesInput'))):
                merged_gtf.append(os.path.join(gtfOutputDir, 'cuffmerge_'+str(i),'merged.gtf'))
            self.setOutput('merged_gtf',merged_gtf)
        else:
            self.setOutput('merged_gtf',None)
        
    def call(self, *args):
        htseqUpstream = args[0]              
        self.setParamIO('assembliesInput',htseqUpstream.getOutput('assembliesOutput'))
        # self.setParamIO('gtfInput1',htseqUpstream.getOutput('gtfOutput1'))
            
    def _singleRun(self, i):
        faInput1 = self.getParamIO('faInput1')
        gtfInput1 = self.getParamIO('gtfInput1')
        assembliesInput = self.getInputList('assembliesInput')
        gtfOutputDir = self.getParamIO('gtfOutputDir')
        cmdline = ['cuffmerge',
                    '-g', gtfInput1,
                    '-s', faInput1,
                    '-o', os.path.join(gtfOutputDir, 'cuffmerge_'+str(i)),
                    '-p', str(self.getParam('threads')),
                    assembliesInput[i]
                    ]
                    
        result = self.callCmdline('V2', cmdline)
        f = open(self.convertToRealPath(os.path.join(Configure.getTmpDir(),'stdout.txt')),'wb')   
        f.write(result.stdout)
        f.close()
            
    def getMarkdownEN(self,):
        mdtext = """
### cuffmerge Result
The cuffmerge result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
Total map reads means that total number of reads mapped to genome
        """.format(mapRs = os.path.join(Configure.getTmpDir(),'stdout.txt'))

        return mdtext
            
        