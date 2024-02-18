'''
Run by:  python3 test_paradis.py paradis_default
'''

import numpy as np
from ctypes import *
import time, sys, os

sys.path.extend([os.path.abspath('../../python'),os.path.abspath('../../lib')])
sys.path.extend([os.path.abspath('../../extensions/paradis/python'),os.path.abspath('../../extensions/paradis/lib')])

home_lib = __import__('Home')
from paradis_util import *

paradis = ParaDiS(home_lib)

'''
Main Program Starts Here
def main():
'''
if True:
    global home, param

    home = paradis.paradis_init("paradis_default")
    param = home.contents.param

    param.contents.timestepIntegrator = b'forceBsubcycle'
    param.contents.subInteg0Integ1 = b'RKF-RKF'

    node_data = np.array([3, 
         0,0,      -4000.0000,   500.0000,        6000.0000,   1,   7,
               0,1,       0.5773503000,     0.5773503000,    -0.5773503000,
                          0.,               1.,               1.,
         0,1,       17.0000,     500.0000,        6000.0000,  2,   0,
               0,0,     -0.5773503000,    -0.5773503000,      0.5773503000,
                        0.,              1.,              1.,
               0,2,      0.5773503000,     0.5773503000,     -0.5773503000,
                        0.,              1.,              1.,
         0,2,       4000.0000,    500.0000,        6000.0000,   1,   7,
               0,1,     -0.5773503000,    -0.5773503000,     0.5773503000,
                        0.,              1.,              1.])

    for iter in range(3):
        print('***************************************************')
        #paradis.RecycleAllNodes(home)
        paradis.FreeAllNodes(home)
        paradis.AddNodesFromArray(home, node_data.ctypes.data_as(POINTER(c_double)))

        t_begin = time.time()  

        maxstep = 10

        for tstep in range(maxstep):
            t0 = time.time()  
            home.contents.cycle = tstep
            paradis.ParadisStep(home)
            t1 = time.time()
            print('step = %d/%d  time used = %g s'%(tstep, maxstep, t1-t0))

        t_end = time.time()  
        print('Total time used = %g s'%(t_end-t_begin))
