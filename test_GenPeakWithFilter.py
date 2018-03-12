# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 11:28
@Author  : Weizhang
@FileName: test_GenPeakWithFilter.py
"""

from StepBase import Configure,Schedule
from GenPeakWithFilter import GenPeakWithFilter

a=GenPeakWithFilter(summitInput='./minidata/atac/BedForTest/output_summits.bed',
                    blacklist='./minidata/atac/BedForTest/consensusBlacklist.bed',
                    topPeak=50000)

Schedule.run()

