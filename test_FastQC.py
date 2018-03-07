
# coding: utf-8
from stepbase import Step,Configure,Schedule
from FastQC import FastQC

# Folder Test
# fastqc = FastQC('./minidata/test_fastqc/','fastq',)
# FileTest
fastqc = FastQC('./minidata/test_fastqc/SRR1294845_2.fastq','fastq',)

Schedule.run()
