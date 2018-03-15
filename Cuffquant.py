# -*- coding: utf-8 -*-
"""
Created on 2018-3-10 09:24:38
@author: Song Shaoming
"""

from StepBase import Step,Configure
import os

class Cuffquant(Step):
	def __init__(self,
				 bamInput = None,
				 gtfInput = None,
				 outputDir = None,

				 threads = None,
				 cmdParam = None,
				 **kwargs
				 ):
		super(Step, self).__init__(cmdParam,**kwargs)

		self.setParamIO('bamInput',bamInput)
		self.setParamIO('gtfInput',gtfInput)
		self.setParamIO('outputDir',outputDir)

		self.initIO()

		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)

		self._setUpstreamSize(2)

	def impInitIO(self,):
		bamInput = self.getParamIO('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		outputDir = self.getParamIO('outputDir')
		if outputDir is None:
			self.setParamIO('outputDir',Configure.getTmpDir())

		if gtfInput is not None:	
			self.setInput('gtfInput',gtfInput)

		self.setInputDirOrFile('bamInput',bamInput)

		#self.setOutputDir1To1('outputDir',outputDir,'cuffquant','suffix','bamInput')
		#self.setOutput('assembliesOutput',os.path.join(Configure.getTmpDir(), 'assembleOfAbundances.txt'))

		if bamInput is not None:
			self._setInputSize(len(self.getInputList('bamInput')))
			abundances_cxb = list()
			for i in range(len(self.getInputList('bamInput'))):
				abundances_cxb.append(os.path.join(outputDir, 'cuffquant_' + str(i), 'abundances.cxb'))
			self.setOutput('abundances_cxb',abundances_cxb)
		else:
			self.setOutput('abundances_cxb',None)

	def call(self, *args):
		bamUpstream = args[0]
		gtfUpstream = args[1]

		self.setParamIO('bamInput',bamUpstream.getOutput('bamOutput'))
		self.setParamIO('gtfInput',gtfUpstream.getOutput('merged_gtf'))

	def _singleRun(self,i):
		bamInput = self.getInputList('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		outputDir = self.getParamIO('outputDir')
		print(outputDir)

		cmdline = [
				'cuffquant',
				'-o',os.path.join(outputDir,'cuffquant_' + str(i)),
				'-p',str(self.getParam('threads')),
				gtfInput[0],
				bamInput[i]
				]
		self.callCmdline('V1',cmdline)