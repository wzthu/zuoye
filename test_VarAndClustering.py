# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 17:41
@Author  : Weizhang
@FileName: test_VarAndClustering.py
"""

from VarAndClustering import VarAndClustering
from StepBase import Configure,Schedule

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')

a=VarAndClustering(bamInput='./minidata/atac/bam_sorted_rmdup',
                   peakInput='./minidata/atac/BedForTest/top_peaks.bed',
                   threads=4,
                   genome='hg19')

Schedule.run()