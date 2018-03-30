from SortBam import SortBam
from StepBase import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

sb = SortBam(bamInput='./step_00_StarAlign/starAligned.out.sam', sortOrder = 'queryname')
Schedule.run()
