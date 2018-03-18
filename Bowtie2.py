# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 16:28:52 2018

@author: WeiZheng
"""

from StepBase import Step,Configure
import subprocess
import os

class Bowtie2(Step):
    """
    Bowtie2 is a Step to align ATAC-seq reads to reference genome.
    See __init__ to initialize this step.
    > Bowtie2(): __init__  parameters
        fastqInput1: str or str list
            mate 1 fastq file path(s) or a directory contain all of fastq files 
        fastqInput2: str or str list
            mate 2 fastq file path(s) or a directory contain all of fastq files
        bt2Idx: str
            bowtie 2 index prefix, for example: /home/data/hg19
        samOutputDir: str
            the output directory of samfiles
        mapRsOutputDir:str
            the mapping result storage directory 
        threads: int
            the threads number will be used
        isdovetail: bool
            If the mates “dovetail”, that is if one mate alignment extends past 
            the beginning of the other such that the wrong mate begins upstream, 
            consider that to be concordant. 
        isNoUnal: bool    
            By default, Bowtie2 do not looks for discordant alignments 
            if it cannot find any concordant alignments.
        isNoMixed: bool
            By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, 
            it then do not tries to find alignments for the individual mates. This option disables that behavior.
        X: int
            The maximum fragment length for valid paired-end alignments.
        cmdParam: str or list of string
            current unsupported    
    > Bowtie2()(): __call__ parameters
    Available upstream objects combinations:
        (AdapterRemoval)
        (SRAToFastq)
        #如果有两个输入，写成(类名1，类名2)，以此类推
    
    """
    def __init__(self,
                 fastqInput1 = None, 
                 fastqInput2 = None, 
                 bt2Idx = None,                 
                 samOutputDir = None, 
                 mapRsOutputDir = None,
                 threads = None,
                 isdovetail = True,
                 isNoDiscordant = True,
                 isNoUnal = True,
                 isNoMixed = True,
                 X = 2000,
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        """
        called by 'Bowtie()'
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function. 
        > Parameters
        fastqInput1: str or str list
            mate 1 fastq file path(s) or a directory contain all of fastq files 
        fastqInput2: str or str list
            mate 2 fastq file path(s) or a directory contain all of fastq files
        bt2Idx: str
            bowtie 2 index prefix, for example: /home/data/hg19
        samOutputDir: str
            the output directory of samfiles
        mapRsOutputDir:str
            the mapping result storage directory 
        threads: int
            the threads number will be used
        isdovetail: bool
            If the mates “dovetail”, that is if one mate alignment extends past 
            the beginning of the other such that the wrong mate begins upstream, 
            consider that to be concordant. 
        isNoUnal: bool    
            By default, Bowtie2 do not looks for discordant alignments 
            if it cannot find any concordant alignments.
        isNoMixed: bool
            By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, 
            it then do not tries to find alignments for the individual mates. This option disables that behavior.
        X: int
            The maximum fragment length for valid paired-end alignments.
        cmdParam: str or list of string
            current unsupported    
        """
        
        # set all input and output parameters
        self.setParamIO('fastqInput1',fastqInput1)
        self.setParamIO('fastqInput2',fastqInput2)
        self.setParamIO('bt2Idx',bt2Idx)         
        self.setParamIO('samOutputDir',samOutputDir)
        self.setParamIO('mapRsOutputDir',mapRsOutputDir)    

        # call self.initIO()
        self.initIO()
            
        #set other parameters
        self.setParam('isNoDiscordant', isNoDiscordant)
        self.setParam('isNoUnal', isNoUnal)
        self.setParam('isNoMixed', isNoMixed)
        self.setParam('isdovetail', isdovetail)
        self.setParam('X', X)
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        
        # obtain all input and output parameters        
        fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput2 = self.getParamIO('fastqInput2')
        bt2Idx = self.getParamIO('bt2Idx')        
        samOutputDir = self.getParamIO('samOutputDir')
        mapRsOutputDir = self.getParamIO('mapRsOutputDir')

        #set all input files        
        self.setInputDirOrFile('fastqInput1',fastqInput1)
        self.setInputDirOrFile('fastqInput2',fastqInput2)      
       
        #some special input from __init__ or configure
        if bt2Idx is None:
            self.setInput('bt2IdxFiles', Configure.getConfig('bt2IdxFiles')) 
            self.setParamIO('bt2Idx', Configure.getConfig('bt2Idx'))
        else:
            suffix = ['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
            bt2IdxFiles = [ bt2Idx + s for s in suffix ]
            self.setInput('bt2IdxFiles', bt2IdxFiles)

        # create output file paths and set
        if samOutputDir is None:
            self.setParamIO('samOutputDir',Configure.getTmpDir())
        if mapRsOutputDir is None:
            self.setParamIO('mapRsOutputDir',Configure.getTmpDir())
        self.setOutputDir1To1('samOutput', samOutputDir,None,'sam','fastqInput1') 
        self.setOutputDir1To1('mapRsOutput',mapRsOutputDir,None,'result.txt','fastqInput1')
        
        # set how many sample are there
        if fastqInput1 is not None:
            self._setInputSize(len(self.getInputList('fastqInput1')))
        
    def call(self, *args):
        """
        called by Bowtie2()(upstreamObj)

        """
        # the first object
        fastqUpstream = args[0]      
        
        # set all required input parameters from upstream object
        self.setParamIO('fastqInput1',fastqUpstream.getOutput('fastqOutput1'))
        self.setParamIO('fastqInput2',fastqUpstream.getOutput('fastqOutput2'))
        #self.setParamIO('bt2Idx',bt2Idx)         
        #self.setParamIO('samOutputDir',samOutputDir)
        #self.setParamIO('mapRsOutputDir',mapRsOutputDir) 

 
            
    def _singleRun(self, i):
        """
        create and execute the command line        
        i is the No. of the sample
        """
        #get all input and output
        fastqInput1 = self.getInputList('fastqInput1')
        fastqInput2 = self.getInputList('fastqInput2')
        samOutput = self.getOutputList('samOutput')
        mapRsOutput = self.getOutputList('mapRsOutput')
        bt2IdxFiles = self.getInput('bt2IdxFiles')            
        #combine the command line
        cmdline = [#'/root/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2',
                'bowtie2',
                '-p', str(self.getParam('threads')),
                self.getBoolParamCmd('--dovetail ', 'isdovetail'),
                self.getBoolParamCmd('--no-discordant ','isNoDiscordant'),
                self.getBoolParamCmd('--no-unal ','isNoUnal'),
                self.getBoolParamCmd('--no-mixed ','isNoMixed'),
                '-X' , str(self.getParam('X')),
                self.getUnsetParams(),
                ' -x %s -q -1 %s -q -2 %s -S %s '%(
                    self.getParamIO('bt2Idx'),
                    #bt2IdxFiles[0],
                    fastqInput1[i],
                    fastqInput2[i],
                    samOutput[i]),
                ]
        
        #run commandline           
        result = self.callCmdline('V1',cmdline,stdoutToLog = True)
        #result = self.callCmdline(cmdline,stdoutToLog = False)
        
        #optional
        f = open(self.convertToRealPath(mapRsOutput[i]),'wb')   
        f.write(result.stdout)
        f.write(result.stderr)
        
    def getMarkdownEN(self,):
        a = """
## Bowtie2 mapping result
Bowtie2 mapping result is shown below:
```{{r, echo=FALSE}}
f<-file("{mapRs}")
readLines(f)
```
        """.format(mapRs=self.getOutput("mapRsOutput")[0])
        return a
        
