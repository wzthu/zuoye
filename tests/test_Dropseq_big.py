from FastqToBam import FastqToBam
from BamMerge import BamMerge
from TagBarcode import TagBarcode
from FilterBam import FilterBam
from TrimAdapter import TrimAdapter
from TrimPolyA import TrimPolyA
from BamToFastq import BamToFastq
from StarAlign import StarAlign
from SortBam import SortBam
from MergeBamAlign import MergeBamAlign
from TagGene import TagGene
from DetectError import DetectError
from DigitalExpression import DigitalExpression
from EasyTreat import EasyTreat
from MonocleQC import MonocleQC
from Monocle_dimreduce_cluster import Monocle_dimreduce_cluster
from StepBase import Configure,Schedule
import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

f2b = FastqToBam(fastqInput1='../data/nuclei_900_fastqs/read1', fastqInput2='../data/nuclei_900_fastqs/read2')
                 #fastqInput1 = './minidata/dropseq/read1', fastqInput2 = './minidata/dropseq/read2')
                 #fastqInput1 = '../data/hgmm_100/read1', fastqInput2 = '../data/hgmm_100/read2')
bm = BamMerge()(f2b)
tbc = TagBarcode(baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)(bm)
tbm = TagBarcode(baseStart = 17, baseEnd = 26, baseQuality = 10,
                barcodeRead = 1, discardRead = True, tagName = 'XM', numBaseBelowQuality = 1)(tbc)
fb = FilterBam(tagReject = 'XQ')(tbm)
ta = TrimAdapter(adapterSeq = 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
                  misMatches = 0, numBases = 5)(fb)
tp = TrimPolyA(misMatches = 0, numBases = 6)(ta)
b2f = BamToFastq()(tp)
sa = StarAlign(outFileNamePrefix='star', genomeDir = '../ref/refdata-cellranger-mm10-2.1.0/star/', threads=16)(b2f)
sb = SortBam(sortOrder = 'queryname')(sa)
mba = MergeBamAlign(refInputDir='../ref/refdata-cellranger-mm10-2.1.0/fasta/',
                secondAlign=False, pairedRun=False)(tp, sb)
tg = TagGene(gtfInput = '../ref/refdata-cellranger-mm10-2.1.0/genes/genes.gtf', tag='GE')(mba)
de = DetectError(numCells=900, primerSeqence='AAGCAGTGGTATCAACGCAGAGTACATGGG')(tg)
dge = DigitalExpression(numCells=900)(de)
#et = EasyTreat()(dge)
mq = MonocleQC()(dge)
mdc = Monocle_dimreduce_cluster(cluster_num=5)(mq) 
Schedule.run()
