from MapDNA import Bowtie2
from RemoveAdapter import AdapterRemoval
from stepbase import Configure,Schedule

Configure.setRefDir('/home/zwei/ref')
Configure.setGenome('hg19')


adrm = AdapterRemoval(fastqInput1='./minidata/atac/end1',fastqInput2='./minidata/atac/end2')
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


