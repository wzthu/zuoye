from cellranger import Cellranger
from stepbase import Configure,Schedule

Configure.setRefDir('/home/zwei/ref')
Configure.setGenome('hg19')

test = Cellranger(fastqInput = '/home/cfeng/data/fastqs/', outputdir='test_cellranger', refile = '/home/cfeng/data/refdata-cellranger-hg19_and_mm10-1.2.0', expectcells=100)
Schedule.run()

print('')
