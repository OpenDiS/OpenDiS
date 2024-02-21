import numpy as np
from ctypes import c_double, POINTER, byref
real8 = c_double

try:
    pydis_lib = __import__('pydis_lib')
    found_pydis = True
except ImportError:
    found_pydis = False
    raise

class paradis_lib:
    def __init__(self):
        for item in pydis_lib.__dict__:
            self.__dict__[item] = pydis_lib.__dict__[item]

    def paradis_init(self):
        """
        Initialize ParaDiS home structure
        """
        home = POINTER(pydis_lib.Home_t)()
        self.ParadisInit_lean(byref(home))

        return home
