# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from ..core import Step,Configure

import os

class SC3_DE(Step):
    def __init__(self,
                 sce_file = None,
                 outputpath = None,
                 cluster_num = 0,
                 cmdParam=None,
                 **kwargs):
        """
        SC3_DE is XXX. Need to use this 
        Step as the downstream of SingleCellExperiment.
        >SC3_DE():_init_parameters
            sce: str
            The R workspace saved by SingleCellExperiment Step.
            outputpath: str
            A str indicates the name of appointed folder that saves outputs.You should
            build that folder in advance. The absolute path is also legel. 
            cluster_num: int
            The number of clusters.
            set to 0 will auto estimate the cluster number
            cmdParam: str or list of string
            current unsupported
        >SC3_DE()():_call_parameters
            Avaliabel upstream objects combinations:
            (SingleCellExperiment)
        """
        super(Step, self).__init__(cmdParam,**kwargs)
        
        

        # set all input and output parameters
        self.setInput('sce_file',sce_file)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters


        self.setParam('cluster_num',cluster_num)
        #self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
               
        outputpath = self.getParamIO('outputpath')  
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath') 

        sce_file = self.getParamIO('sce_file') 
        # create output file paths and set
        self.setOutputDir1To1('outputpaths',None,outputpath, 'Expression','jpg','sce_file')
        self.setOutputDir1To1('fastqcOutput_Expreesion_jpg',None,outputpath, 'Expression.','jpg','sce_file')
        self.setOutputDir1To1('fastqcOutput_Expreesion_jpg',None,outputpath, 'Expression.','jpg','sce_file')

       
    def call(self,*args):

        Upstream = args[0]
        if isinstance(Upstream,FastqDump):
            fastqInput = Upstream.getOutput('fastqOutput1')
            fastqInput.extend( Upstream.getOutput('fastqOutput2'))
            self.setParamIO('fastqInput', fastqInput)
        # set all required input parameters from upstream object
        #上游可能為 “fastqInput1”,“fastqInput2”,“fastqOnput1”
        #self.setParamIO('fastqInput',Upstream.getOutput('fastqOutput1'))

        print("Call the UpStream Node, Not implementation")

        #other things

    def getMarkdownEN(self,):
        mdtext = """
## FastQC Usage

FastQC('/path/to/input_fastq',[fastq/sra],'/path/to/output_dir')  
Attention!  
* In this release, only .fastq format file can be setting as input!  

## FastQC Quality Control Result  
The FastQC Quality Control is shown below:  
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(knitr)
library(kableExtra)
fastq_name <- {fastq_name}
fastqc_html <- {fastqc_report}
fq <- cbind(fastq_name, fastqc_html)
colnames(fq) <- c("Fastq Name", "Fastq Report")
kable(fq, "html") %>% kable_styling() %>% scroll_box(width = "1100px", height = "500px")
```

        
"""
 
        return mdtext
            
            
    def _singleRun(self,i):
        # obtain all input and output dir list
        fastqInput = self.getInputList('fastqInput')
        fastqcOutputDir = self.getParamIO('fastqcOutputDir')

        cmdline =['fastqc',
                   '-t',str(self.getParam('threads')),
                   '-o',fastqcOutputDir,
                   fastqInput[i]
                   ]
        self.callCmdline('V1', cmdline)