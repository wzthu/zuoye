from MergeBamAlign import MergeBamAlign
from StepBase import Configure,Schedule

import os

mba = MergeBamAlign(unmappedBamInput='./minidata/dropseq/tmp/unaligned_mc_tagged_polyA_filtered.bam', alignedBamInput='./minidata/dropseq/tmp/aligned.sorted.bam',
                bamOutputDir='./minidata/dropseq/tmp', refSequence='../../dropseq/refdata-cellranger-hg19_and_mm10-2.1.0/fasta/genome.fa',
                secondAlign=False, pairedRun=False)
Schedule.run()
