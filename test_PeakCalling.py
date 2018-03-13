# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/12 21:52
@Author  : Weizhang
@FileName: test_PeakCalling.py
"""

from PeakCalling import PeakCalling
from StepBase import Configure, Schedule

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')


test=PeakCalling(bedInput='./minidata/atac/others/MergedAllFiles.bed')


Schedule.run()


