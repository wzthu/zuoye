# -*- coding: utf-8 -*-
"""
Created on 2018-3-8 20:06:30

@author: Song Shaoming
"""

from stepbase import Step,Configure
import os

class Bamsort(Step):
	def __init__(self,
				 bamInput = None,
				 bamOutputDir = None,
				 cmdParam = None,
				 **kwargs
				 ):
		super(Step, self).__init__(cmdParam,**kwargs)

		self.setParamIO('bamInput',bamInput)
		self.setParamIO('bamOutputDir',bamOutputDir)

		self.initIO()

	def impInitIO(self,):
		bamInput = self.getParamIO('bamInput')
		bamOutputDir = self.getParamIO('bamOutputDir')

		self.setInputDirOrFile('bamInput',bamInput)
		self.setOutputDir1To1('bamOutput',bamOutputDir,'sorted','bam','bamInput')

		if bamInput is not None:
			self._setInputSize(len(self.getInputList('bamInput')))


	def call(self, *args):

		bamUpstream = args[0]
		self.setParamIO('bamInput',bamUpstream.getOutput('bamOutput'))

	def _singleRun(self, i):
		bamInput = self.getInputList('bamInput')
		bamOutput = self.getOutputList('bamOutput')

		cmdline = [
				'samtools sort',
				bamInput[i],
				'-o',bamOutput[i]
				]
		self.callCmdline(cmdline)
