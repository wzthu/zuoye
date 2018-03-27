from StepBase import Schedule,Configure
from FastqDump import FastqDump
from Hisat2 import Hisat2
from SamToBam import SamToBam
from BamSort import BamSort
from Cufflinks import Cufflinks
from Cuffmerge import Cuffmerge
from Cuffquant import Cuffquant
from Cuffdiff import Cuffdiff

Configure.setIdentity('sqchen0327')
def smartseq_flow(sraInput, ht2Idx_ref, gtf_ref, fa_ref, threads):

    fastq_dump = FastqDump(sraInput1=sraInput)

    hisat = Hisat2(ht2Idx=ht2Idx_ref)(fastq_dump)

    sam2bam = SamToBam(threads=threads)(hisat)

    bamsort = BamSort()(sam2bam)

    cufflinks = Cufflinks(gtfInput=gtf_ref,threads=threads)(bamsort)

    cuffmerge = Cuffmerge(faInput1=fa_ref,gtfInput1=gtf_ref,threads=threads)(cufflinks)

    cuffquant = Cuffquant(threads=threads)(bamsort,cuffmerge)

    cuffdiff = Cuffdiff(faInput=fa_ref,threads=threads)(cuffmerge,cuffquant)

    Schedule.run()

def smartseq_flow2(sraInput, ht2Idx_ref, gtf_ref, fa_ref, threads):

    fastq_dump = FastqDump(sraInput1=sraInput)

    hisat = Hisat2(ht2Idx=ht2Idx_ref)(fastq_dump)

    sam2bam = SamToBam(threads=threads)(hisat)

    bamsort = BamSort()(sam2bam)

    cufflinks = Cufflinks(gtfInput=gtf_ref,threads=threads)(bamsort)

    cuffmerge = Cuffmerge(faInput1=fa_ref,gtfInput1=gtf_ref,threads=threads)(cufflinks)

    cuffquant = Cuffquant(threads=threads)(bamsort,cuffmerge)

    cuffnorm = Cuffnorm(threads=threads)(cuffquant,cuffmerge)

    Schedule.run()

smartseq_flow(sraInput='/data8t_1/chenshengquan/minidata/test_sra', 
    ht2Idx_ref="/data8t_1/ref/smartseq/hg19_index/genome", 
    gtf_ref='/data8t_1/ref/smartseq/genome.gtf', 
    fa_ref='/data8t_1/ref/smartseq/hg19.fa', 
    threads=16)
smartseq_flow2(sraInput='/data8t_1/chenshengquan/minidata/test_sra', 
    ht2Idx_ref="/data8t_1/ref/smartseq/hg19_index/genome", 
    gtf_ref='/data8t_1/ref/smartseq/genome.gtf', 
    fa_ref='/data8t_1/ref/smartseq/hg19.fa', 
    threads=16)