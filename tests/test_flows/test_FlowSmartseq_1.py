from hcacn.core import Configure,Schedule
from hcacn.flows import FlowSmartseq_1, FlowSmartseq_1

Configure.setIdentity('sqchen_f1')
rs=FlowSmartseq_1(sraInput='/data8t_1/chenshengquan/minidata/test_sra',refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultDir='./resultFlowSmartseq_1')()
# Configure.setIdentity('sqchen_f2')
# rs=FlowSmartseq_1(sraInput='/data8t_1/chenshengquan/minidata/test_sra',refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultDir='./resultFlowSmartseq_2')()
