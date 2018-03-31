# -*- coding: utf-8 -*-

from hcacn.core import Configure, Schedule
from hcacn.steps import SingleCellExperiment


#Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
#Configure.setGenome('hg19')
Configure.setIdentity('yinqijin')

stf = SingleCellExperiment(matrix_file='/data8t_1/hca/zuoye/minidata/downstream/matrix/matrix.csv',
			ann_file = '/data8t_1/hca/zuoye/minidata/downstream/annotation/annotation.csv',
                 matrix_format = 'ORIGIN',
                 #outputpath = None,
			)



Schedule.run()
