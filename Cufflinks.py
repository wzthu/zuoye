# -*- coding: utf-8 -*-
"""
Created on 2018-3-7 10:44:26

@author: Song Shaoming
"""

from StepBase import Step,Configure
import os

class Cufflinks(Step):
	def __init__(self,
				 bamInput = None,
				 gtfInput = None,
				# fragBiasCorrectInput = None,
				 outputDir = None,
				 
				 threads = None,
				 ismultiReadCorrect = None,
				 isupperQuartileForm = None,
				 istotalHitsNorm = True,
				 fragLenMean = 200,
				 fragLenStdDev = 80,
				 cmdParam = None,
				 **kwargs
				):
		super(Step, self).__init__(cmdParam,**kwargs)

		self.setParamIO('bamInput',bamInput)
		self.setParamIO('gtfInput',gtfInput)
		#self.setParamIO('fragBiasCorrectInput',fragBiasCorrectInput)
		self.setParamIO('outputDir',outputDir)
		self.initIO()

		self.setParam('ismultiReadCorrect',ismultiReadCorrect)
		self.setParam('fragLenMean',fragLenMean)
		self.setParam('fragLenStdDev',fragLenStdDev)
		self.setParam('isupperQuartileForm',isupperQuartileForm)
		self.setParam('istotalHitsNorm',istotalHitsNorm)

		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)

	def impInitIO(self,):
		bamInput = self.getParamIO('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		outputDir = self.getParamIO('outputDir')
		#fragBiasCorrectInput = self.getParamIO('fragBiasCorrectInput')

		self.setInputDirOrFile('bamInput',bamInput)
		self.setInputDirOrFile('gtfInput',gtfInput)
		#self.setInputDirOrFile('fragBiasCorrectInput',fragBiasCorrectInput)

		#if fragBiasCorrectInput is None:
		#	self.setParamIO('fragBiasCorrectInput',' ')

		self.setOutputDirNTo1('outputDir',outputDir,'genes.fpkm_tracking','bamInput')

		if bamInput is not None:
			self._setInputSize(len(self.getInputList('bamInput')))

	def call(self, *args):

		bamUpstream = args[0]

		self.setParamIO('bamInput',bamUpstream.getOutput('bamOutput'))

	def _singleRun(self,i):
		bamInput = self.getInputList('bamInput')
		gtfInput = self.getInputList('gtfInput')
		#fragBiasCorrectInput = self.getInputList('fragBiasCorrectInput')
		print('uhhh')
		outputDir = self.getOutputList('outputDir')

		cmdline = [
				'cufflinks',
				'-p',str(self.getParam('threads')),
				self.getBoolParamCmd('-u','ismultiReadCorrect'),
				self.getBoolParamCmd('-N','isupperQuartileForm'),
				self.getBoolParamCmd('--total-hits-norm','istotalHitsNorm'),
				'-m',str(self.getParam('fragLenMean')),
				'-s',str(self.getParam('fragLenStdDev')),
				'-G',gtfInput[i],
				'-o',outputDir[i],
				bamInput[i]
				]
		self.callCmdline(cmdline)

