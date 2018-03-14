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
		#self.setOutputDir1To1('outputDir',outputDir,'Cuffnorm','suffix','cxbInput')
		if cxbInput is not None:
			self._setInputSize(len(self.getInputList('cxbInput')))
			isoforms_fpkm_tracking=list()
			genes_fpkm_tracking=list()
			cds_fpkm_tracking=list()
			tss_groups_fpkm_tracking=list()	
			isoforms_count_tracking=list()
			genes_count_tracking=list()
			cds_count_tracking=list()
			tss_groups_count_tracking=list()
			for i in range(len(self.getInputList('cxbInput'))):
				isoforms_fpkm_tracking.append(os.path.join(outputDir,'cuffnorm_'+str(i),'isoforms_fpkm_tracking'))
				genes_fpkm_tracking.append(os.path.join(outputDir, 'cuffnorm_'+str(i),'genes.fpkm_tracking'))
				cds_fpkm_tracking.append(os.path.join(outputDir, 'cuffnorm_'+str(i),'cds.fpkm_tracking'))
				tss_fpkm_tracking.append(os.path.join(outputDir, 'cuffnorm_'+str(i),'tss_groups.fpkm_tracking'))
				isoforms_count_tracking.append(os.path.join(outputDir,'cuffnorm_'+str(i),'isoforms_count_tracking'))
				genes_count_tracking.append(os.path.join(outputDir, 'cuffnorm_'+str(i),'genes.count_tracking'))
				cds_count_tracking.append(os.path.join(outputDir, 'cuffnorm_'+str(i),'cds.count_tracking'))
				tss_count_tracking.append(os.path.join(outputDir, 'cuffnorm_'+str(i),'tss_groups.count_tracking'))
			self.setOutput('isoforms_fpkm_tracking',isoforms_fpkm_tracking)
			self.setOutput('genes_fpkm_tracking',genes_fpkm_tracking)
			self.setOutput('cds_fpkm_tracking',cds_fpkm_tracking)
			self.setOutput('tss_groups_fpkm_tracking',tss_groups_fpkm_tracking)
			self.setOutput('isoforms_count_tracking',isoforms_count_tracking)
			self.setOutput('genes_count_tracking',genes_count_tracking)
			self.setOutput('cds_count_tracking',cds_count_tracking)
			self.setOutput('tss_groups_count_tracking',tss_groups_count_tracking)
		else:
			self.setOutput('isoforms_fpkm_tracking',None)
			self.setOutput('genes_fpkm_tracking',None)
			self.setOutput('cds_fpkm_tracking',None)
			self.setOutput('tss_groups_fpkm_tracking',None)
			self.setOutput('isoforms_count_tracking',None)
			self.setOutput('genes_count_tracking',None)
			self.setOutput('cds_count_tracking',None)
			self.setOutput('tss_groups_count_tracking',None)

	  def call(self, *args):
		cxbUpstream = args[0]
		fd = open(cxbUpstream.getOutput('assembliesOutput'),'r')
		cxb = ''
		marker = ''
		for lines in fd.readlines():
			content = lines.split('/')
			cxb = cxb + lines.strip() + ','
			marker = marker + content[-2] + ','
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
				'-o',os.path.join(outputDir,'cuffnorm_'+str(i),cxbInput[i]),
				'-p',str(self.getParam('threads')),
				'-L',markerInput,
				gtfInput,
				cxbInput[i]
				]
		self.callCmdline('V1',cmdline,stdoutToLog = True)
