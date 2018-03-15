# -*- coding: utf-8 -*-
"""
Created on 2018-3-12 16:48:36

@author: Song Shaoming
"""

from StepBase import Step,Configure
import os
import subprocess

class Cuffdiff(Step):
	def __init__(self,
				 faInput = None,
				 gtfInput = None,
				 cxbInput = None,
				 markerInput = None,
				 outputDir = None,

				 threads = None,
				 cmdParam = None,
				 **kwargs
				):
		super(Step, self).__init__(cmdParam,**kwargs)

		self.setParamIO('faInput',faInput)
		self.setParamIO('gtfInput',gtfInput)
		self.setParamIO('cxbInput',cxbInput)
		self.setParamIO('markerInput',markerInput)
		self.setParamIO('outputDir',outputDir)
		self.initIO()
		
		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)
		self._setUpstreamSize(2)


	def impInitIO(self,):
		faInput = self.getParamIO('faInput')
		gtfInput = self.getParamIO('gtfInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		if cxbInput is not None:
			for i,item in enumerate(cxbInput.split(' ')):
				self.setInput('cxbInput_%d'%i,item)

		if gtfInput is not None:
			self.setInput('gtfInput',gtfInput)

		if faInput is None:
			faInput = Configure.getConfig('')
			self.setInput('faInput',faInput)
			self.setParamIO('faInput',faInput)
		else:
			self.setInput('faInput',faInput)
			
		if outputDir is None:
			self.setParamIO('outputDir',Configure.getTmpDir())

		outputDir = self.getParamIO('outputDir')

		isoforms_fpkm_table=list()
		genes_fpkm_table=list()
		cds_fpkm_table=list()
		tss_groups_fpkm_table=list()	
		isoforms_count_table=list()
		genes_count_table=list()
		cds_count_table=list()
		tss_groups_count_table=list()
		isoforms_exp_diff=list()
		genes_exp_diff=list()
		cds_exp_diff=list()
		tss_groups_exp_diff=list()

		isoforms_fpkm_table.append(os.path.join(outputDir,'cuffdiff_'+str(0),'isoforms.fpkm_table'))
		genes_fpkm_table.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'genes.fpkm_table'))
		cds_fpkm_table.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'cds.fpkm_table'))
		tss_groups_fpkm_table.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'tss_groups.fpkm_table'))
		isoforms_count_table.append(os.path.join(outputDir,'cuffdiff_'+str(0),'isoforms.count_table'))
		genes_count_table.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'genes.count_table'))
		cds_count_table.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'cds.count_table'))
		tss_groups_count_table.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'tss_groups.count_table'))
		isoforms_exp_diff.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'isoforms_exp.diff'))
		genes_exp_diff.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'genes_exp.diff'))
		cds_exp_diff.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'cds_exp.diff'))
		tss_groups_exp_diff.append(os.path.join(outputDir, 'cuffdiff_'+str(0),'tss_groups_exp.diff'))

		self.setOutput('isoforms_fpkm_table',isoforms_fpkm_table)
		self.setOutput('genes_fpkm_table',genes_fpkm_table)
		self.setOutput('cds_fpkm_table',cds_fpkm_table)
		self.setOutput('tss_groups_fpkm_table',tss_groups_fpkm_table)
		self.setOutput('isoforms_count_table',isoforms_count_table)
		self.setOutput('genes_count_table',genes_count_table)
		self.setOutput('cds_count_table',cds_count_table)
		self.setOutput('tss_groups_count_table',tss_groups_count_table)
		self.setOutput('isoforms_exp_diff',isoforms_exp_diff)
		self.setOutput('genes_exp_diff',genes_exp_diff)
		self.setOutput('cds_exp_diff',cds_exp_diff)
		self.setOutput('tss_groups_exp_diff',tss_groups_exp_diff)

		self._setInputSize(1)

	def call(self, *args):
		cxbUpstream = args[1]
		gtfUpstream = args[0]

		cxb = cxbUpstream.getOutput('abundances_cxb')
		#print('===================')
		#print(cxb)		
		marker = ''
		for i in range(len(cxb)):
			#cxb[i] = self.convertToRealPath(cxb[i])
			marker = marker + cxb[i].strip().split('/')[-2] + ','
		cxb = ' '.join(cxb)
		marker = marker[:-1]
		self.setParamIO('cxbInput',cxb)
		self.setParamIO('markerInput',marker)

		self.setParamIO('gtfInput',gtfUpstream.getOutput('merged_gtf'))
		print('==================================')
		print(gtfUpstream.getOutput('merged_gtf'))

	def _singleRun(self,i):
		gtfInput = self.getParamIO('gtfInput')
		faInput = self.getParamIO('faInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		cxbFinPath = []
		for i,item in enumerate(cxbInput.split(' ')):
			cxbFinPath.append(self.getInput('cxbInput_'+str(i)))
		cxbFin = ' '.join(cxbFinPath)

		cmdline = [
				'cuffdiff',
				'-o',os.path.join(outputDir,'cuffdiff_' + str(0)),
				'-p',str(self.getParam('threads')),
				'-L',markerInput,
				'-b',faInput,
				
				gtfInput[0],
				cxbFin
				]
		self.callCmdline('V1',cmdline,stdoutToLog = True)	