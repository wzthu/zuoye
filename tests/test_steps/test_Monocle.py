# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
 """
from hcacn.core import Configure,Schedule
from hcacn.steps import MonocleQC, MonocleDC


Configure.setIdentity('zywang')
#Configure.enableDocker(False)


#Configure.enableDocker(False)
MonocleQC_result = MonocleQC(matrixdata='out_gene_exon_tagged.dge.txt', outputpath='outputresult') 

# To see if all input and output parameter are right in paramsIO 
MonocleQC_result.paramsIO

# To see if other parameters are right in params
MonocleQC_result.params

# To see if all input files are right
MonocleQC_result.inputs 

MonocleQC_result.outputs

MonocleDC_result = MonocleDC(outputpath='outputresult')(MonocleQC_result)

Schedule.run()


