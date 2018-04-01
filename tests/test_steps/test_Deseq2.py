# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
 """
from hcacn.core import Configure, Schedule
from hcacn.steps import Deseq2

Configure.setIdentity('zywang')
Deseq2(matrixdata = "Deseq2testdata.txt", annotation = "condition.csv", outputpath = None)
Schedule.run()