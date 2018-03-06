from Mapping import Bowtie2
from RemoveAdapter import AdapterRemoval
from stepbase import Configure,Schedule

Configure.setRefDir('/home/zwei/ref')
Configure.setGenome('hg19')


adrm = AdapterRemoval(fastqInput1=['./minidata/atac/sample1/chr20_1.1.fq',
                                   './minidata/atac/sample2/chr20_1.2.fq'],
                      fastqInput2=['./minidata/atac/sample1/chr20_2.1.fq',
                                   './minidata/atac/sample2/chr20_2.2.fq'])
"""
# To see if all input and output parameter are right in paramsIO 
adrm.paramsIO

# To see if other parameters are right in params
adrm.params

# To see if all input files are right
adrm.inputs 

# To see if all input files are right
adrm.params
"""

rs=Bowtie2()(adrm)

"""
bt2Obj = Bowtie2()
rs=Bowtie2()(adrm)
"""


Schedule.run()


