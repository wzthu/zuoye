# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 17:41
@Author  : Weizhang
@FileName: test_VarAndClustering.py
"""

from VarAndClustering import VarAndClustering
from StepBase import Configure,Schedule

Configure.setRefDir('/home/hca/zhangwei1/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

a=VarAndClustering(bamInput='/home/hca/zhangwei1/bam_sorted_rmdup',
                   peakInput='./step_12_GenPeakWithFilter/MergedAllFiles_summits_filterd.bed',
                   threads=4,
                   genome='hg19')

Schedule.run()