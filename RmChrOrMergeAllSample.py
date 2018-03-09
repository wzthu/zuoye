# -*- coding: utf-8 -*-
"""

designed for scATAC-seq

@author: Weizhang

Remove some chromatin(optional) or merge all sample bed files to a big bed file.

"""

from StepBase import Step, Configure


class RmChrOrMergeAllSample(Step):
    def __init__(self,
                 bedInput=None,
                 bedOutputDir=None,
                 savedchr = None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bedInput', bedInput)
        self.setParamIO('bedOutputDir', bedOutputDir)

        self.initIO()
        self.setParam('savedchr', savedchr)


    def impInitIO(self):
        bedInput = self.getParamIO('bedInput')
        bedOutputDir = self.getParamIO('bedOutputDir')

        # set all input files
        self.setInputDirOrFile('bedInput', bedInput)
        # set all output files
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, 'bed', 'bedInput')
        self.setOutputDirNTo1('mergedfilename', None, 'MergedAllFiles.bed', 'bedInput')
        if bedInput is not None:
            self._setInputSize(len(self.getInputList('bedInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bedInput', samUpstream.getOutput('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bedInput = self.getInputList('bedInput')
        bedOutput = self.getOutputList('bedOutput')
        mergedfilename = self.getOutputList('mergedfilename')[0]
        chr_info = self.getParam('savedchr')

        if chr_info is None:  # only merge, without remove chromatin
            g = open(r'%s' % mergedfilename, 'w')  # for saving merged reads
            for i in range(len(bedInput)):
                f = open(r'%s' % bedInput[i], 'r')  # for reading every sample
                for line in f.readlines():
                    temp = line.split()
                    if temp[0] in chr_info:
                        new_line = '\t'.join([str(i) for i in temp])
                        g.write(new_line + '\n')

                f.close()
            g.close()
        else:  # merge and remove chromatin(saving file for every sample)
            g = open(r'%s' % mergedfilename, 'w')  # for saving merged reads
            for i in range(len(bedInput)):
                f = open(r'%s' % bedInput[i], 'r')  # for reading every sample
                g1 = open(r'%s' % bedOutput[i], 'w')  # for saving sample reads
                for line in f.readlines():
                    temp = line.split()
                    if temp[0] in chr_info:
                        new_line = '\t'.join([str(i) for i in temp])
                        g.write(new_line + '\n')
                        g1.write(new_line + '\n')
                g1.close()
                f.close()
            g.close()




