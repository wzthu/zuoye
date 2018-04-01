# -*- coding: utf-8 -*-

from hcacn.core import Configure, Schedule
from hcacn.steps import SingleCellExperiment,SC3_DE,SC3_Cluster


#Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
#Configure.setGenome('hg19')
Configure.setIdentity('yinqijin')


# test single input file
# stf = SingleCellExperiment(matrix_file='/data8t_1/hca/zuoye/minidata/downstream/matrix/matrix.csv',
# 			ann_file = '/data8t_1/hca/zuoye/minidata/downstream/annotation/annotation.csv',
#                  matrix_format = 'ORIGIN',
#                  #outputpath = None,
# 			)

# test multi input file and setting output folder
sce = SingleCellExperiment(matrix_file='/data8t_1/hca/zuoye/minidata/downstream/matrix/',
			ann_file = '/data8t_1/hca/zuoye/minidata/downstream/annotation/',
                 matrix_format = 'ORIGIN',
                 outputpath = 'step_my',
			)

sc3_de = SC3_DE(
 	outputpath = None,
    cluster_num = 4,
                 )(sce)


sc3_cluster = SC3_Cluster(
	outputpath = "./step_my_sc3_cluster",
    cluster_num = 0,
                 )(sce)
Schedule.run()
