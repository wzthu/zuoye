from Cellranger import Cellranger
from Seurat import Seurat
from StepBase import Configure, Schedule

test = Cellranger(fastqInput = '/home/cfeng/data/test/',  refile = '/home/cfeng/data/refdata-cellranger-hg19_and_mm10-1.2.0', expectcells=100)
test2 = Seurat(rscript = '/home/cfeng/test_celranger/Seurat.R')(test)
Schedule.run()

print('')

