# -*- coding: utf-8 -*-

from Bowtie2 import Bowtie2
from AdapterRemoval import AdapterRemoval
from StepBase import Configure,Schedule
from SamToBam import SamToBam
from BamSort import BamSort
from RmDuplicates import RmDuplicates
from BamToBed import BamToBed
from RmChrOrMergeAllSample import RmChrOrMergeAllSample
from MergeToFrag import MergeToFrag
from BedSort import BedSort
from SRAToFastq import SRAToFastq

Configure.setRefDir('/home/hca/zhangwei/hg19')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

stf = SRAToFastq(sraInput='./minidata/atac/SraForTest')

print("SRAToFastq output:")
print(stf.getOutput('fastqOutput1'))
print(stf.getOutput('fastqOutput2'))
print("\n")

adrm = AdapterRemoval()(stf)

print("AdapterRemoval output:")
print(adrm.getOutput('fastqOutput1'))
print(adrm.getOutput('fastqOutput2'))
print("\n")

rs=Bowtie2()(adrm)

print("Bowtie2 output:")
print(rs.getOutput('samOutput'))
print("\n")


sb=SamToBam(threads=6)(rs)

print("SamToBam output:")
print(sb.getOutput('bamOutput'))
print("\n")

bs=BamSort(threads=3)(sb)

print("BamSort output:")
print(bs.getOutput('bamOutput'))
print("\n")

rd=RmDuplicates(picard='/home/hca/zhangwei/software1/picard.jar', memory='-Xmx4g')(bs)

print("RmDuplicates output:")
print(bs.getOutput('bamOutput'))
print("\n")


# bb=BamToBed()(rd)
#
# print("BamToBed output:")
# print(bb.getOutput('bedOutput'))
# print("\n")
#
# # just for a test
# mtf1=MergeToFrag()(bb)
#
# print("MergeToFrag1 output:")
# print(mtf1.getOutput('bedOutput'))
# print("\n")
#
# chr_info = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
#             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
#
# rcmas=RmChrOrMergeAllSample(savedchr=chr_info)(bb)
#
# print("RmChrOrMergeAllSample output:")
# print(rcmas.getOutput('mergedfilename'))
# print(rcmas.getOutput('bedOutput'))
# print("\n")
#
# mtf2=MergeToFrag()(rcmas)
#
# print("MergeToFrag2 output:")
# print(mtf2.getOutput('bedOutput'))
# print("\n")
#
# bedsort=BedSort()(rcmas)
#
# print("BedSort output:")
# print(bedsort.getOutput('bedOutput'))
# print("\n")
#

Schedule.run()








#
# from Bowtie2 import Bowtie2
# from AdapterRemoval import AdapterRemoval
# from StepBase import Configure,Schedule
# from SamToBam import SamToBam
# from BamSort import BamSort
# from RmDuplicates import RmDuplicates
# from BamToBed import BamToBed
# from RmChrOrMergeAllSample import RmChrOrMergeAllSample
# from MergeToFrag import MergeToFrag
# from BedSort import BedSort
# from SRAToFastq import SRAToFastq
#
# Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
# Configure.setGenome('hg19')
#
# stf=SRAToFastq(sraInput='./minidata/atac/SraForTest')
#
# print("SRAToFastq output:")
# print(stf.getOutput('fastqOutput1'))
# print(stf.getOutput('fastqOutput2'))
# print("\n")
#
# adrm = AdapterRemoval()(stf)
#
# print("AdapterRemoval output:")
# print(adrm.getOutput('fastqOutput1'))
# print(adrm.getOutput('fastqOutput2'))
# print("\n")
#
# rs=Bowtie2()(adrm)
#
# print("Bowtie2 output:")
# print(rs.getOutput('samOutput'))
# print("\n")
#
# sb=SamToBam(threads=6)(rs)
#
# print("SamToBam output:")
# print(sb.getOutput('bamOutput'))
# print("\n")
#
# bs=BamSort(threads=3)(sb)
#
# print("BamSort output:")
# print(bs.getOutput('bamOutput'))
# print("\n")
#
# rd=RmDuplicates(picard='/home/wzhang/software/Picard/picard.jar', memory='-Xmx4g')(bs)
#
# print("RmDuplicates output:")
# print(bs.getOutput('bamOutput'))
# print("\n")
#
# bb=BamToBed()(rd)
#
# print("BamToBed output:")
# print(bb.getOutput('bedOutput'))
# print("\n")
#
# # just for a test
# mtf1=MergeToFrag()(bb)
#
# print("MergeToFrag1 output:")
# print(mtf1.getOutput('bedOutput'))
# print("\n")
#
# chr_info = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
#             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
#
# rcmas=RmChrOrMergeAllSample(savedchr=chr_info)(bb)
#
# print("RmChrOrMergeAllSample output:")
# print(rcmas.getOutput('mergedfilename'))
# print(rcmas.getOutput('bedOutput'))
# print("\n")
#
# mtf2=MergeToFrag()(rcmas)
#
# print("MergeToFrag2 output:")
# print(mtf2.getOutput('bedOutput'))
# print("\n")
#
# bedsort=BedSort()(rcmas)
#
# print("BedSort output:")
# print(bedsort.getOutput('bedOutput'))
# print("\n")
#
# Schedule.run()





