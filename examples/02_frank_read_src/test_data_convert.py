import numpy as np
import sys, os

pyexadis_paths = ['../../python', '../../lib', '../../core/exadis/python','../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from pyexadis_base import ExaDisNet
from pydis import DisNet

import pyexadis
pyexadis.initialize()

# run the simulation defined in test_frank_read_src_exadis.py
import test_frank_read_src_exadis
test_frank_read_src_exadis.main()

# access the global variables defined in test_frank_read_src_exadis.py
net = test_frank_read_src_exadis.net
sim = test_frank_read_src_exadis.sim

G0 = net.get_disnet(ExaDisNet)
G1 = net.get_disnet(DisNet)

# sanity check
G1_sanity_check = G1.is_sane()
if G1_sanity_check:
    print("1. test G1 sanity check" + '\033[32m' + " PASSED" + '\033[0m')
else:
    print("1. test G1 sanity check" + '\033[31m' + " FAILED" + '\033[0m')

# test export_data and import data 
G2 = DisNet()
G2.import_data(G1.export_data())

G1_eq_G2 = G1.is_equivalent(G2) and G2.is_equivalent(G1)
if G1_eq_G2:
    print("2. test G1_eq_G2" + '\033[32m' + " PASSED" + '\033[0m')
else:
    print("2. test G1_eq_G2" + '\033[31m' + " FAILED" + '\033[0m')


G3 = G1.to_networkx()
G4 = DisNet()
G4.from_networkx(G3)

# test conversion to and from networkx
G1_eq_G4 = G1.is_equivalent(G4) and G4.is_equivalent(G1)
if G1_eq_G4:
    print("3. test G1_eq_G4" + '\033[32m' + " PASSED" + '\033[0m')
else:
    print("3. test G1_eq_G4" + '\033[31m' + " FAILED" + '\033[0m')
