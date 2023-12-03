'''
Run by:  python3 test_paradis.py paradis_default
'''

from ctypes import *
import time, sys, os

sys.path.extend([os.path.abspath('../../extensions/paradis/python'),os.path.abspath('../../extensions/paradis/lib')])

use_GPU = True if len(sys.argv)> 2 and sys.argv[2]=='use_GPU' else False

home_lib = __import__('Home_gpu' if use_GPU else 'Home')
from paradis_util import *

paradis = ParaDiS(home_lib)

'''
Main Program Starts Here
'''
def main():
    global home, param, use_GPU

    taskname = sys.argv[1]
    home = paradis.paradis_init(taskname)
    param = home.contents.param

    print("use_GPU = ", use_GPU)

    param.contents.timestepIntegrator = b'forceBsubcycle'
    if use_GPU:
        param.contents.subInteg0Integ1 = b'GPU'
    else:
        param.contents.subInteg0Integ1 = b'RKF-RKF'

    if use_GPU: paradis.InitializeParadisGPU(home)

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


if __name__ == "__main__":
     main()
