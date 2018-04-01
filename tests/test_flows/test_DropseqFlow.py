# -*- coding: utf-8 -*-
from hcacn.core import Configure,Schedule
from hcacn.flows import DropseqFlow
Configure.setIdentity('lcy')

df=DropseqFlow(fastqInput1='./minidata/dropseq/read1', fastqInput2='./minidata/dropseq/read2',refDir=None, refdir='/data8t_1/ref/dropseq/', genome='hg19-and-mm10')()

