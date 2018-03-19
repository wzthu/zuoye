
# coding: utf-8
from StepBase import Step,Configure,Schedule
from FastQC import FastQC
Configure.setIdentity("yinqijin")

Configure.enableDocker(True)
# Folder Test
# fastqc = FastQC('./minidata/test_fastqc/','fastq',)
# FileTest
fastqc = FastQC('./minidata/smartseq/fastq/','fastq',)

Schedule.run()
