import networkx as nx
import numpy as np

from dd_proto import DDSim

ddsim = DDSim()

def set_param():
    print("set_param")
    ddsim.param.mu = 160e9
    ddsim.param.nu = 0.31

def init_loop(radius=1.0, burg=None, plane_normal=None):
    print("init_loop")
    # create loop in ddsim.disnet 
    return None

def print_force():
    print("print_force")
    return None

def main():
    set_param()
    init_loop()

    print("calculate node forces")
    ddsim.NodeForce()

    print_force()


if __name__ == "__main__":
    main()
