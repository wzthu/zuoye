# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 11:28
@Author  : Weizhang
@FileName: test_GenPeakWithFilter.py
"""

from StepBase import Configure,Schedule
from GenPeakWithFilter import GenPeakWithFilter

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

a=GenPeakWithFilter(summitInput='./minidata/atac/others/output_summits.bed',
                    blacklist='./minidata/atac/others/consensusBlacklist.bed',
                    topPeak=50000)

Schedule.run()

