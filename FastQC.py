
# coding: utf-8
from StepBase import Step,Configure,Schedule
import os

class FastQC(Step):

    def  __init__(self,

                  fastqInput=None,
                  fileFormat =None,
                  fastqcOutputDir = None,
                  threads = None,
                  cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)

        # set all input and output parameters
        self.setParamIO('fastqInput',fastqInput)
        if fastqcOutputDir == None:
            self.setParamIO('fastqcOutputDir',Configure.getTmpDir())
        else:
            self.setParamIO('fastqcOutputDir',fastqcOutputDir)

        # call self.initIO()
        self.initIO()

        #set other parameters
        #self.setParam('isNoDiscordant', isNoDiscordant)
        self.setParam('fileFormat',fileFormat)
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)

        print (self.params)

    def impInitIO(self,):

        # obtain all input and output parameters
        fastqInput = self.getParamIO('fastqInput')
        fastqcOutputDir = self.getParamIO('fastqcOutputDir')

        #set all input files
        self.setInputDirOrFile('fastqInput',fastqInput)


        # create output file paths and set
        self.setOutputDir1To1('fastqcOutput', fastqcOutputDir,'fastqcTest','_fastqc','fastqInput')



        # set how many sample are there
        if fastqInput is not None:
            self._setInputSize(len(self.getInputList('fastqInput')))


    def call(self,*args):

        Upstream = args[0]

        # set all required input parameters from upstream object
        #上游可能為 “fastqInput1”,“fastqInput2”,“fastqOnput1”
        #self.setParamIO('fastqInput',Upstream.getOutput('fastqOutput1'))

        print("Call the UpStream Node, Not implementation")

        #other things

    def _singleRun(self,i):
        # obtain all input and output dir list
        fastqInput = self.getInputList('fastqInput')
        fastqcOutput = self.getOutputList('fastqcOutput')



        if not os.path.exists(fastqcOutput[i]):
            os.mkdir(fastqcOutput[i])


        cmdline = ['fastqc',
                   '-t',str(self.getParam('threads')),
                   '-o',fastqcOutput[i],
                   fastqInput[i]
                  ]
        self.callCmdline('V1', cmdline)
