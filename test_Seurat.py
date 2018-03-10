from StepBase import Configure,Schedule
from Seurat import Seurat
import os

test = Seurat(outputdir = 'test_seurat', rscript = '/data/test_celranger/Seurat.R')
Schedule.run()

print('')
