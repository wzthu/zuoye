from stepbase import Configure,Schedule
from Seurat import Seurat

Configure.setRefDir('/home/zwei/ref')



adrm = AdapterRemoval(fastqInput1='./minidata/atac/sample1/chr20_1.1.fq',fastqInput2='./minidata/atac/sample1/chr20_2.1.fq')
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

rs=Bowtie2()
results = rs(adrm)

"""
bt2Obj = Bowtie2()
rs=Bowtie2()(adrm)
"""
sb = SamToBam()(rs)

Schedule.run()


