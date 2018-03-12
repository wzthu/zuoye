# -*- coding: utf-8 -*-
"""
Created on 2018-3-10 15:05:34

@author: Song Shaoming
"""

from stepbase import Step,Configure
import os

class Cuffnorm(Step):
	def __init__(self,
				 gtfInput = None,
				 cxbInput = None,
				 outputDir = None,
				 markerInput = None,

				 threads = None,
				 cmdParam = None,
				 **kwargs
				 ):
		super(Step, self).__init__(cmdParam,**kwargs)

		self.setParamIO('gtfInput',gtfInput)
		self.setParamIO('cxbInput',cxbInput)
		self.setParamIO('markerInput',markerInput)
		self.setParamIO('outputDir',outputDir)
		self.initIO()

		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)

	def impInitIO(self,):
		gtfInput = self.getParamIO('gtfInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		self.setInputDirOrFile('cxbInput',cxbInput)
		self.setOutputDir1To1('outputDir',outputDir,'Cuffnorm','suffix','cxbInput')
		if cxbInput is not None:
			self._setInputSize(len(self.getInputList('cxbInput')))

	def call(self, *args):
		cxbUpstream = args[0]
		fd = open(cxbUpstream.getOutput('assembliesOutput'),'r')
		cxb = ''
		marker = ''
		for lines in fd.readlines():
			content = lines.split('/')
			size = len(content)
			tag = '/'
			for i in range(size - 1):
				tag = tag + content[i] + '/'
			tag = tag + 'abundunces.cxb'
			cxb = cxb + tag + ','
			marker = marker + content[3] + ','
		cxb = cxb[:-1]
		marker = marker[:-1]
		self.setParamIO('cxbInput',cxb)
		self.setParamIO('markerInput',marker)

	def _singleRun(self,i):
		gtfInput = self.getParamIO('gtfInput')
		cxbInput = self.getInputList('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getOutputList('outputDir')

		cmdline = [
				'Cuffnorm',
				'-o',outputDir[i],
				'-p',str(self.getParam('threads')),
				'-L',markerInput,
				gtfInput,
				cxbInput[i]
				]
		self.callCmdline('V1',cmdline,stdoutToLog = True)