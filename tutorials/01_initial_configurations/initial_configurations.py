import numpy as np
import sys, os

option = sys.argv[1] if len(sys.argv) > 1 else None

if option == "1":
    print('Create simple dislocation source using pydis')
    
    pydis_paths = ['../../python', '../../lib', '../../core/pydis/python']
    [sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

    from framework.disnet_manager import DisNetManager
    from pydis import DisNode, DisNet, Cell, VisualizeNetwork

    Ldis = 100.0 # dislocation length
    Lbox = 2*Ldis # simulation box size
    burg = 1.0/np.sqrt(2.0)*np.array([1.,1.,0.]) # Burgers vector
    plane = np.array([-1.,1.,1.]) # plane normal

    # Simulation cell object
    cell = Cell(h=Lbox*np.eye(3), is_periodic=[True,True,True])

    # List of nodes, nodes = [x,y,z,constraint]
    linevec = Ldis*burg # line vector
    rn = np.array([[*(-0.5*linevec), DisNode.Constraints.PINNED_NODE],
                   [*( 0.0*linevec), DisNode.Constraints.UNCONSTRAINED],
                   [*( 0.5*linevec), DisNode.Constraints.PINNED_NODE]])
    rn[:,0:3] += cell.center() # translate line to the center of the box

    # List of segments, links = [node1,node2,burg,plane]
    links = np.array([[0, 1, *burg, *plane],
                      [1, 2, *burg, *plane]])

    G = DisNet(cell=cell, rn=rn, links=links)
    N = DisNetManager(G)

    VisualizeNetwork().plot_disnet(N, block=True)


elif option == "2":
    print('Create simple dislocation source using pyexadis')
    
    pyexadis_paths = ['../../python', '../../core/exadis/python']
    [sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
    
    from framework.disnet_manager import DisNetManager
    import pyexadis
    from pyexadis_base import ExaDisNet, NodeConstraints, VisualizeNetwork
    pyexadis.initialize()

    Ldis = 100.0 # dislocation length
    Lbox = 2*Ldis # simulation box size
    burg = 1.0/np.sqrt(2.0)*np.array([1.,1.,0.]) # Burgers vector
    plane = np.array([-1.,1.,1.]) # plane normal

    # Simulation cell object
    cell = pyexadis.Cell(h=Lbox*np.eye(3), is_periodic=[True,True,True])

    # List of nodes, nodes = [x,y,z,constraint]
    linevec = Ldis*burg # line vector
    nodes = np.array([[*(-0.5*linevec), NodeConstraints.PINNED_NODE],
                      [*( 0.0*linevec), NodeConstraints.UNCONSTRAINED],
                      [*( 0.5*linevec), NodeConstraints.PINNED_NODE]])
    nodes[:,0:3] += cell.center() # translate line to the center of the box

    # List of segments, segs = [node1,node2,burg,plane]
    segs = np.array([[0, 1, *burg, *plane],
                     [1, 2, *burg, *plane]])

    G = ExaDisNet(cell, nodes, segs)
    N = DisNetManager(G)
    
    VisualizeNetwork().plot_disnet(N, block=True)
    
    if not sys.flags.interactive:
        pyexadis.finalize()


elif option == "3":
    print('Read initial configuration from ParaDiS file')
    
    pyexadis_paths = ['../../python', '../../core/exadis/python']
    [sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
    
    from framework.disnet_manager import DisNetManager
    import pyexadis
    from pyexadis_base import ExaDisNet, VisualizeNetwork
    pyexadis.initialize()
    
    G = ExaDisNet().read_paradis('../../examples/10_strain_hardening/180chains_16.10e.data')
    N = DisNetManager(G)
    
    VisualizeNetwork().plot_disnet(N, block=True)
    
    if not sys.flags.interactive:
        pyexadis.finalize()


elif option == "4":
    print('Create initial configuration of infinite lines')
    
    pyexadis_paths = ['../../python', '../../core/exadis/python']
    [sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
    
    from framework.disnet_manager import DisNetManager
    import pyexadis
    from pyexadis_base import ExaDisNet, VisualizeNetwork
    pyexadis.initialize()
    
    G = ExaDisNet().generate_line_config(
        crystal='fcc', Lbox=1000.0, num_lines=12, theta=90.0, maxseg=100.0, seed=1234
    )
    N = DisNetManager(G)
    
    VisualizeNetwork().plot_disnet(N, block=True)
    
    from pyexadis_utils import write_data, write_vtk
    write_data(N, 'infinite_lines.data')
    write_vtk(N, 'infinite_lines.vtk')
    
    if not sys.flags.interactive:
        pyexadis.finalize()


elif option == "5":
    print('Create initial configuration of prismatic loops')
    
    pyexadis_paths = ['../../python', '../../core/exadis/python']
    [sys.path.append(os.path.abspath(path)) for path in pyexadis_paths if not path in sys.path]
    
    from framework.disnet_manager import DisNetManager
    import pyexadis
    from pyexadis_base import ExaDisNet, VisualizeNetwork
    pyexadis.initialize()
    
    G = ExaDisNet().generate_prismatic_config(
        crystal='bcc', Lbox=1000.0, num_loops=12, radius=100.0, maxseg=50.0, seed=5678, uniform=True
    )
    N = DisNetManager(G)
    
    VisualizeNetwork().plot_disnet(N, block=True)
    
    if not sys.flags.interactive:
        pyexadis.finalize()
    
    
elif option is None:
    print(f'Specify one of the following options:')
    print(' 1: Create simple dislocation source using pydis')
    print(' 2: Create simple dislocation source using pyexadis')
    print(' 3: Read initial configuration from ParaDiS file')
    print(' 4: Create initial configuration of infinite lines')
    print(' 5: Create initial configuration of prismatic loops')
    
else:
    print(f'Unknown option {option}')
