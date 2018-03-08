# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""
from stepbase import Configure,Schedule
from Seurat import Seurat
import os

Seurat_result = Seurat(inputDir='F:/HCA/mm10', outputDir='F:/HCA/stepone') 

# To see if all input and output parameter are right in paramsIO 
Seurat_result.paramsIO

# To see if other parameters are right in params
Seurat_result.params

# To see if all input files are right
Seurat_result.inputs 

Schedule.run()


