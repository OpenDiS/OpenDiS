Python binding: ``pyexadis``
============================

This section documents the raw ExaDiS classes and functions available through the python interface, ``pyexadis``.
These modules are direct bindings to the backend C++ modules implemented in ExaDiS.
For documentation of the ExaDiS python modules available in OpenDiS format, see :doc:`Python modules: pyexadis_base <python_modules>`.
For documentation about the backend C++ modules please see the :doc:`Developer guide <../../developer_guide/index>` section of the documentation.

.. toctree::
   :maxdepth: 1


General Attributes
------------------

.. py:attribute:: pyexadis.BCC_CRYSTAL
   :type: int
   :noindex:

   Constant for BCC crystal type.

.. py:attribute:: pyexadis.FCC_CRYSTAL
   :type: int
   :noindex:

   Constant for FCC crystal type.

General Classes
---------------

.. py:class:: pyexadis.Params

   Simulation parameters.

   .. py:method:: __init__(crystal="", burgmag, mu, nu, a, maxseg, minseg, rann=-1.0, rtol=-1.0, maxdt=1e-7, nextdt=1e-12, split3node=1)
      :noindex:

      Initialize parameters.

      :param crystal: Crystal type (string)
      :param burgmag: Burgers vector magnitude
      :param mu: Shear modulus
      :param nu: Poisson's ratio
      :param a: Dislocation core radius
      :param maxseg: Maximum line discretization length
      :param minseg: Minimum line discretization length
      :param rann: Annihilation distance
      :param rtol: Error tolerance
      :param maxdt: Maximum timestep size
      :param nextdt: Starting timestep size
      :param split3node: Enable splitting of 3-nodes

   .. py:attribute:: crystalparams
      :type: CrystalParams
      :noindex:

      Crystal parameters.

   .. py:attribute:: burgmag
      :type: float
      :noindex:

      Burgers vector magnitude (scaling length).

   .. py:attribute:: mu
      :type: float
      :noindex:

      Shear modulus.

   .. py:attribute:: nu
      :type: float
      :noindex:

      Poisson's ratio.

   .. py:attribute:: a
      :type: float
      :noindex:

      Dislocation core radius.

   .. py:attribute:: maxseg
      :type: float
      :noindex:

      Maximum line discretization length.

   .. py:attribute:: minseg
      :type: float
      :noindex:

      Minimum line discretization length.

   .. py:attribute:: rann
      :type: float
      :noindex:

      Annihilation distance.

   .. py:attribute:: rtol
      :type: float
      :noindex:

      Error tolerance.

   .. py:attribute:: maxdt
      :type: float
      :noindex:

      Maximum timestep size.

   .. py:attribute:: nextdt
      :type: float
      :noindex:

      Starting timestep size.

   .. py:attribute:: split3node
      :type: int
      :noindex:

      Enable splitting of 3-nodes.

---

.. py:class:: pyexadis.CrystalParams

   Parameters for crystal type, orientation, and glide planes.
   
   .. py:method:: set_crystal_type(str)
      :noindex:

      Set the crystal type.

   .. py:attribute:: R
      :type: array
      :noindex:

      Crystal orientation matrix.

   .. py:attribute:: use_glide_planes
      :type: bool
      :noindex:

      Use and maintain dislocation glide planes.

   .. py:attribute:: enforce_glide_planes
      :type: bool
      :noindex:

      Enforce glide planes option.

   .. py:attribute:: num_bcc_plane_families
      :type: int
      :noindex:

      Number of BCC plane families (1, 2, or 3).

---

.. py:class:: pyexadis.Crystal

   Crystal type and orientation.

   .. py:method:: __init__(type)
      :noindex:

      Initialize with crystal type.

   .. py:method:: __init__(type, R)
      :noindex:

      Initialize with crystal type and orientation matrix.

   .. py:method:: __init__(crystalparams)
      :noindex:

      Initialize with crystal parameters.

   .. py:attribute:: type
      :type: int
      :noindex:

      Index of the crystal type.

   .. py:attribute:: R
      :type: array
      :noindex:

      Crystal orientation matrix.

   .. py:method:: set_orientation(R)
      :noindex:

      Set crystal orientation matrix.

   .. py:method:: set_orientation(euler_angles)
      :noindex:

      Set crystal orientation via Euler angles.

---

.. py:class:: pyexadis.Cell

   Simulation cell and periodic boundary conditions.

   .. py:method:: __init__(Lbox, centered=False)
      :noindex:

      Initialize cubic cell.

   .. py:method:: __init__(Lvecbox, centered=False)
      :noindex:

      Initialize cell with vector box.

   .. py:method:: __init__(bmin, bmax)
      :noindex:

      Initialize cell with bounds.

   .. py:method:: __init__(h, origin=Vec3(0.0), is_periodic=[PBC_BOUND, PBC_BOUND, PBC_BOUND])
      :noindex:

      Initialize cell with matrix, origin, and periodicity.

   .. py:method:: __init__(cell)
      :noindex:

      Copy constructor.

   .. py:attribute:: h
      :type: array
      :noindex:

      Cell matrix.

   .. py:attribute:: origin
      :type: array
      :noindex:

      Origin of the cell.

   .. py:method:: center()
      :noindex:

      Returns the center of the cell.

   .. py:method:: closest_image(Rref, R)
      :noindex:

      Returns the closest image of an array of positions from another reference.

   .. py:method:: pbc_fold(R)
      :noindex:

      Fold an array of positions to the primary cell.

   .. py:method:: is_inside(R)
      :noindex:

      Checks if a position is inside the primary cell.

   .. py:method:: are_inside(R)
      :noindex:

      Checks if an array of positions are inside the primary cell.

   .. py:method:: is_triclinic()
      :noindex:

      Returns if the box is triclinic.

   .. py:method:: is_periodic()
      :noindex:

      Get the cell periodic boundary condition flags along the 3 dimensions.

   .. py:method:: get_bounds()
      :noindex:

      Get the (orthorhombic) bounds of the cell.

   .. py:method:: volume()
      :noindex:

      Returns the volume of the cell.

---

.. py:class:: pyexadis.NodeTag

   Tag identifying a node by domain and local index.

   .. py:method:: __init__(domain, index)
      :noindex:

      :param domain: Domain index (int)
      :param index: Local index (int)

   .. py:attribute:: domain
      :type: int
      :noindex:

      Domain index.

   .. py:attribute:: index
      :type: int
      :noindex:

      Local index.

---

.. py:class:: pyexadis.DisNode

   Dislocation node, including position, constraint, and tag.

   .. py:method:: __init__(pos, constraint)
      :noindex:

      :param pos: Node position (Vec3)
      :param constraint: Node constraint flag (int)

   .. py:attribute:: tag
      :type: NodeTag
      :noindex:

      Node tag (domain, index).

   .. py:attribute:: constraint
      :type: int
      :noindex:

      Node constraint flag.

   .. py:attribute:: pos
      :type: Vec3
      :noindex:

      Node position (x, y, z).

---

.. py:class:: pyexadis.DisSeg

   Dislocation segment between two nodes.

   .. py:method:: __init__(n1, n2, burg, plane=Vec3(0.0))
      :noindex:

      :param n1: Segment start node index (int)
      :param n2: Segment end node index (int)
      :param burg: Segment Burgers vector (Vec3)
      :param plane: Segment plane normal (Vec3, default zero vector)

   .. py:attribute:: n1
      :type: int
      :noindex:

      Segment start node index.

   .. py:attribute:: n2
      :type: int
      :noindex:

      Segment end node index.

   .. py:attribute:: burg
      :type: Vec3
      :noindex:

      Segment Burgers vector.

   .. py:attribute:: plane
      :type: Vec3
      :noindex:

      Segment plane normal.

---

.. py:class:: pyexadis.Conn

   Connectivity information for a node.

   .. py:method:: __init__()
      :noindex:

      Default constructor.

   .. py:attribute:: num
      :type: int
      :noindex:

      Number of connections.

   .. py:method:: node(i)
      :noindex:

      Get the node index of the i-th connection.

      :param i: Connection index (int)
      :rtype: int

   .. py:method:: seg(i)
      :noindex:

      Get the segment index of the i-th connection.

      :param i: Connection index (int)
      :rtype: int

   .. py:method:: order(i)
      :noindex:

      Get the order value of the i-th connection.

      :param i: Connection index (int)
      :rtype: int

   .. py:method:: add_connection(node, seg, order)
      :noindex:

      Add a connection.

      :param node: Node index (int)
      :param seg: Segment index (int)
      :param order: Order value (int)
      :rtype: bool

   .. py:method:: remove_connection(i)
      :noindex:

      Remove the i-th connection.

      :param i: Connection index (int)

---

.. py:class:: pyexadis.SerialDisNet

   Serial dislocation network object. This is the direct binding to the C++ implementation. In the python interface, a wrapper ExaDisNet object is usually used instead as a proxy for a SerialDisNet object.

   .. py:method:: __init__(cell)
      :noindex:

      :param cell: Simulation cell (Cell)

   .. py:method:: number_of_nodes()
      :noindex:

      Get the number of nodes.

      :rtype: int

   .. py:method:: number_of_segs()
      :noindex:

      Get the number of segments.

      :rtype: int

   .. py:method:: dislocation_density(burgmag)
      :noindex:

      Get the dislocation density.

      :rtype: float

   .. py:attribute:: cell
      :type: Cell
      :noindex:

      Simulation cell.

   .. py:method:: nodes(i)
      :noindex:

      Get the i-th node by reference.

      :param i: Node index (int)
      :rtype: DisNode

   .. py:method:: segs(i)
      :noindex:

      Get the i-th segment by reference.

      :param i: Segment index (int)
      :rtype: DisSeg

   .. py:method:: conn(i)
      :noindex:

      Get the connectivity object for the i-th node by reference.

      :param i: Node index (int)
      :rtype: Conn

   .. py:method:: find_connection(n1, n2)
      :noindex:

      Find connection index c in the Conn array of node 1 that connects to node 2, i.e. where Conn(n1).node(c) == n2.

      :param n1: Node 1 index (int)
      :param n2: Node 2 index (int)
      :rtype: int

   .. py:method:: generate_connectivity()
      :noindex:

      Generate internal connectivity array for the network.

   .. py:method:: _add_node(pos, constraint=UNCONSTRAINED)
      :noindex:

      Add node (x, y, z[, constraint]) to the network.

      :param pos: Position (Vec3)
      :param constraint: Constraint flag (int, default UNCONSTRAINED)

   .. py:method:: _add_seg(n1, n2, burg, plane=Vec3(0.0))
      :noindex:

      Add segment (n1, n2, burg[, plane]) to the network.

      :param n1: Start node index (int)
      :param n2: End node index (int)
      :param burg: Burgers vector (Vec3)
      :param plane: Plane normal (Vec3, default zero vector)

   .. py:method:: move_node(node_index, new_pos, dEp)
      :noindex:

      Move a node to a new position.

      :param node_index: Index of node (int)
      :param new_pos: New position (Vec3)
      :param dEp: Reference to the plastic strain variable to be incremented (Mat33)

   .. py:method:: split_seg(seg_index, new_pos)
      :noindex:

      Split a segment by inserting a new node at new_pos.

      :param seg_index: Segment index (int)
      :param new_pos: New node position (Vec3)
      :returns: Index of the new node (int)

   .. py:method:: split_node(node_index, arms)
      :noindex:

      Split a node given the list of arms to be transferred to the new node.

      :param node_index: Node index (int)
      :param arms: Node arms (array)
      :returns: Index of the new node (int)

   .. py:method:: merge_nodes(node_index1, node_index2, dEp)
      :noindex:

      Merge two nodes into the first node.

      :param node_index1: First node index (int)
      :param node_index2: Second node index (int)
      :param dEp: Reference to the plastic strain variable to be incremented (Mat33)
      :returns: Merge error (bool)

   .. py:method:: merge_nodes_position(node_index1, node_index2, new_pos, dEp)
      :noindex:

      Merge two nodes into the first node at the new position.

      :param node_index1: First node index (int)
      :param node_index2: Second node index (int)
      :param new_pos: New position (Vec3)
      :param dEp: Reference to the plastic strain variable to be incremented (Mat33)
      :returns: Merge error (bool)

   .. py:method:: remove_segs(seg_indices)
      :noindex:

      Remove segments by indices.

      :param seg_indices: List of segment indices

   .. py:method:: remove_nodes(node_indices)
      :noindex:

      Remove nodes by indices.

      :param node_indices: List of node indices

   .. py:method:: purge_network()
      :noindex:

      Purge the network (remove unused nodes/segs).

   .. py:method:: update()
      :noindex:

      Update network memory after modifications.

---

.. py:class:: pyexadis.ExaDisNet

   Dislocation network object, supporting import/export, manipulation, and analysis of network data.

   .. py:method:: __init__()
      :noindex:

      Default constructor.

   .. py:method:: __init__(cell, nodes, segs)
      :noindex:

      Construct network from cell, nodes, and segments.

      :param cell: Cell object
      :param nodes: List of node data (list of lists: [dom, id, x, y, z, constraint])
      :param segs: List of segment data (list of lists: [n1, n2, burg, plane])

   .. py:method:: import_data(cell, nodes, segs)
      :noindex:

      Set the network with (cell, nodes, segs) data.

      :param cell: Cell object
      :param nodes: List of node data
      :param segs: List of segment data

   .. py:method:: number_of_nodes()
      :noindex:

      Returns the number of nodes in the network.

      :rtype: int

   .. py:method:: number_of_segs()
      :noindex:

      Returns the number of segments in the network.

      :rtype: int

   .. py:method:: is_sane()
      :noindex:

      Checks if the network is sane (valid connectivity, no corrupt data).

      :rtype: bool

   .. py:method:: get_cell()
      :noindex:

      Get the cell containing the network by value.

      :rtype: Cell

   .. py:method:: get_nodes_array()
      :noindex:

      Get the list of nodes in the network.

      :returns: List of [dom, id, x, y, z, constraint] for each node

   .. py:method:: get_segs_array()
      :noindex:

      Get the list of segments in the network.

      :returns: List of [n1, n2, burg, plane] for each segment

   .. py:method:: get_forces()
      :noindex:

      Get the list of node forces (fx, fy, fz) in the network.

      :returns: List of [fx, fy, fz] for each node

   .. py:method:: get_velocities()
      :noindex:

      Get the list of node velocities (vx, vy, vz) in the network.

      :returns: List of [vx, vy, vz] for each node

   .. py:method:: set_positions(positions)
      :noindex:

      Set the list of node positions (x, y, z) in the network.

      :param positions: List of [x, y, z] for each node

   .. py:method:: set_forces(forces)
      :noindex:

      Set the list of node forces (fx, fy, fz) in the network.

      :param forces: List of [fx, fy, fz] for each node

   .. py:method:: set_velocities(velocities)
      :noindex:

      Set the list of node velocities (vx, vy, vz) in the network.

      :param velocities: List of [vx, vy, vz] for each node

   .. py:method:: write_data(filename)
      :noindex:

      Write network data in ParaDiS format.

      :param filename: Output file path (str)

   .. py:method:: get_plastic_strain()
      :noindex:

      Returns the plastic strain, plastic spin, and dislocation density, as computed since the last integration step.

      :returns: Plastic strain array, plastic spin array, dislocation density

   .. py:method:: physical_links()
      :noindex:

      Returns the list of segments for each physical dislocation link.

      :returns: List of segment indices for each link

   .. py:method:: _get_crystal()
      :noindex:

      Get the Crystal object by reference.

      :rtype: Crystal

   .. py:method:: _get_serial_network()
      :noindex:

      Get the underlying SerialDisNet object by reference.

      :rtype: SerialDisNet

---

.. py:class:: pyexadis.System(pyexadis.ExaDisNet)

   ExaDiS simulation system, combining a network and simulation parameters.

   .. py:method:: __init__(net, params)
      :noindex:

      Construct a system from a network and parameters.

      :param net: ExaDisNet object
      :param params: Params object

   .. py:method:: set_neighbor_cutoff(cutoff)
      :noindex:

      Set neighbor cutoff of the system.

      :param cutoff: Cutoff value (float)

   .. py:method:: set_applied_stress(stress)
      :noindex:

      Set applied stress of the system.

      :param stress: List or array of [xx, yy, zz, yz, xz, xy]

   .. py:method:: print_timers(dev=False)
      :noindex:

      Print simulation timers.

      :param dev: If True, print device timers (bool, default False)

---

General Functions
-----------------

.. py:function:: pyexadis.initialize(num_threads=-1, device_id=0)

   Initialize the python binding module.

   :param num_threads: Number of OMP threads (default: -1)
   :param device_id: Device GPU index (default: 0)

.. py:function:: pyexadis.finalize()

   Finalize the python binding module.

.. py:function:: pyexadis.read_paradis(filename)

   Read ParaDiS data file.

   :param filename: Path to ParaDiS data file
   :rtype: ExaDisNet

.. py:function:: pyexadis.generate_prismatic_config(crystal, Lbox, numsources, radius, maxseg=-1, seed=1234, uniform=False)
   :noindex:

   Generate a configuration made of prismatic loops in a cubic cell.

   :param crystal: Crystal object
   :param Lbox: Box size
   :param numsources: Number of sources
   :param radius: Loop radius
   :param maxseg: Maximum segment length (default: -1)
   :param seed: Random seed (default: 1234)
   :param uniform: Uniform distribution (default: False)
   :rtype: ExaDisNet

.. py:function:: pyexadis.generate_prismatic_config(crystal, cell, numsources, radius, maxseg=-1, seed=1234, uniform=False)

   Generate a configuration made of prismatic loops in a custom cell.

   :param crystal: Crystal object
   :param cell: Cell object
   :param numsources: Number of sources
   :param radius: Loop radius
   :param maxseg: Maximum segment length (default: -1)
   :param seed: Random seed (default: 1234)
   :param uniform: Uniform distribution (default: False)
   :rtype: ExaDisNet


Force
-----

Force Modules and Factories
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:class:: pyexadis.Force.Force

   ExaDiS base force class.

   .. py:method:: name()
      :noindex:

      Get force name.

---

.. py:class:: pyexadis.Force.CORE_SELF_PKEXT(pyexadis.Force.Force)

   Core, self, and external PK force model.

   .. py:class:: Params
      :noindex:

      Parameters for CORE_SELF_PKEXT model.

      .. py:method:: __init__(Ecore=-1.0, Ecore_junc_fact=1.0)
         :noindex:

         :param Ecore: Core energy (float, default -1.0)
         :param Ecore_junc_fact: Junction energy factor (float, default 1.0)

---

.. py:class:: pyexadis.Force.ForceFFT(pyexadis.Force.Force)

   FFT-based force model.

   .. py:class:: Params
      :noindex:

      Parameters for ForceFFT model.

      .. py:method:: __init__(Ngrid)
         :noindex:

         :param Ngrid: Grid dimensions (list of 3 ints)

   .. py:method:: make(params, fparams, cell)
      :noindex:

      Factory to create a ForceFFT model.

      :param params: General simulation parameters
      :param fparams: ForceFFT.Params object
      :param cell: Cell object
      :rtype: ForceBind

   .. py:method:: get_neighbor_cutoff()
      :noindex:

      Returns the neighbor cutoff used in ForceFFT.

   .. py:method:: get_rcgrid()
      :noindex:

      Returns the grid cutoff used in ForceFFT.

   .. py:method:: interpolate_stress(R)
      :noindex:

      Interpolates stress at query positions from FFT grid values.

      :param R: Array of positions

   .. py:method:: export_stress_gridval()
      :noindex:

      Exports grid stress values.

---

.. py:class:: pyexadis.Force.LINE_TENSION_MODEL(pyexadis.Force.Force)

   Line tension force model.

   .. py:method:: make(params, coreparams)
      :noindex:

      Factory to create a LINE_TENSION_MODEL.

      :param params: Global simulation parameters
      :param coreparams: CORE_SELF_PKEXT.Params object
      :rtype: ForceBind

---

.. py:class:: pyexadis.Force.CUTOFF_MODEL(pyexadis.Force.Force)

   Cutoff-based force model.

   .. py:class:: Params
      :noindex:

      Parameters for CUTOFF_MODEL.

      .. py:method:: __init__(coreparams, cutoff)
         :noindex:

         :param coreparams: CORE_SELF_PKEXT.Params object
         :param cutoff: Cutoff distance (float)

   .. py:method:: make(params, fparams)
      :noindex:

      Factory to create a CUTOFF_MODEL.

      :param params: Global simulation parameters
      :param fparams: CUTOFF_MODEL.Params object
      :rtype: ForceBind

---

.. py:class:: pyexadis.Force.DDD_FFT_MODEL(pyexadis.Force.Force)

   DDD-FFT force model.

   .. py:class:: Params
      :noindex:

      Parameters for DDD_FFT_MODEL.

      .. py:method:: __init__(coreparams, Ngrid)
         :noindex:

         :param coreparams: CORE_SELF_PKEXT.Params object
         :param Ngrid: Grid dimensions (list of 3 ints)

   .. py:method:: make(params, fparams, cell)
      :noindex:

      Factory to create a DDD_FFT_MODEL.

      :param params: Global simulation parameters
      :param fparams: DDD_FFT_MODEL.Params object
      :param cell: Cell object
      :rtype: ForceBind

---

.. py:class:: pyexadis.Force.SUBCYCLING_MODEL(pyexadis.Force.Force)

   Subcycling force model.

   .. py:class:: Params
      :noindex:

      Parameters for SUBCYCLING_MODEL.

      .. py:method:: __init__(coreparams, Ngrid, drift=False, flong_group0=True)
         :noindex:

         :param coreparams: CORE_SELF_PKEXT.Params object
         :param Ngrid: Grid dimensions (list of 3 ints)
         :param drift: Enable drift subcycling scheme (bool, default False)
         :param flong_group0: Integrate long-range force in group 0 (bool, default True)

   .. py:method:: make(params, fparams, cell)
      :noindex:

      Factory to create a SUBCYCLING_MODEL.

      :param params: Global simulation parameters
      :param fparams: SUBCYCLING_MODEL.Params object
      :param cell: Cell object
      :rtype: ForceBind

---

.. py:class:: pyexadis.Force.ForcePython(pyexadis.Force.Force)

   Python-based force model.

   .. py:method:: make(params, force)
      :noindex:

      Factory to create a Python-based force model.

      :param params: Global simulation parameters
      :param force: CalForce python object implementing force methods
      :rtype: ForceBind

---

Force Binder
^^^^^^^^^^^^

.. py:class:: pyexadis.ForceBind

   Wrapper for force computation and related operations.

   .. py:attribute:: force
      :type: Force
      :noindex:

      Internal force object.

   .. py:attribute:: neighbor_cutoff
      :type: float
      :noindex:

      Neighbor cutoff distance.

   .. py:method:: _pre_compute(system)
      :noindex:

      Pre-compute force of the system (internal use).

   .. py:method:: _compute(system)
      :noindex:

      Compute force of the system (internal use).

   .. py:method:: _get_force_fft()
      :noindex:

      Return internal ForceFFT object if it exists in the internal force model.

      :rtype: ForceFFT

   .. py:method:: compute_force(net, applied_stress=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], pre_compute=True)
      :noindex:

      Wrapper to compute nodal forces of the system.

      :param net: ExaDisNet object
      :param applied_stress: Applied stress tensor (list of 6 floats)
      :param pre_compute: Whether to run pre-computation (bool)

   .. py:method:: pre_compute_force(net)
      :noindex:

      Wrapper to perform pre-computations before compute_node_force().

      :param net: ExaDisNet object

   .. py:method:: compute_node_force(net, i, applied_stress)
      :noindex:

      Wrapper to compute the force on a single node.

      :param net: ExaDisNet object
      :param i: Node index (int)
      :param applied_stress: Applied stress tensor (list of 6 floats)

---

Force Calculation Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: pyexadis.compute_force_n2(net, mu, nu, a)

   Compute elastic forces using brute-force N^2 calculation.

   :param net: ExaDisNet object
   :param mu: Shear modulus (float)
   :param nu: Poisson ratio (float)
   :param a: Core cutoff (float)

.. py:function:: pyexadis.compute_force_cutoff(net, mu, nu, a, cutoff, maxseg=0.0)

   Compute elastic forces using a segment pair cutoff.

   :param net: ExaDisNet object
   :param mu: Shear modulus (float)
   :param nu: Poisson ratio (float)
   :param a: Core cutoff (float)
   :param cutoff: Cutoff distance (float)
   :param maxseg: Max segment length (float, default 0.0)

.. py:function:: pyexadis.compute_force_segseglist(net, mu, nu, a, segseglist)

   Compute elastic forces given a list of segment pairs.

   :param net: ExaDisNet object
   :param mu: Shear modulus (float)
   :param nu: Poisson ratio (float)
   :param a: Core cutoff (float)
   :param segseglist: List of segment pairs

---

Mobility
--------

Mobility Binder
^^^^^^^^^^^^^^^

.. py:class:: pyexadis.Mobility

   Mobility law wrapper.

   .. py:method:: compute(system)
      :noindex:

      Compute mobility of the system.

---

Mobility Law Factories
^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: pyexadis.make_mobility_glide(params, mobparams)

   Instantiate a GLIDE mobility law.

.. py:function:: pyexadis.make_mobility_bcc_0b(params, mobparams)

   Instantiate a BCC_0B mobility law.

.. py:function:: pyexadis.make_mobility_bcc_nl(params, mobparams)

   Instantiate a BCC_NL mobility law.

.. py:function:: pyexadis.make_mobility_fcc_0(params, mobparams)

   Instantiate a FCC_0 mobility law.

.. py:function:: pyexadis.make_mobility_fcc_0_fric(params, mobparams)

   Instantiate a FCC_0_FRIC mobility law.

.. py:function:: pyexadis.make_mobility_fcc_0b(params, mobparams)

   Instantiate a FCC_0B mobility law.

.. py:function:: pyexadis.make_mobility_python(params, mobility)

   Instantiate a python-based mobility model.

---

Mobility Calculation Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: pyexadis.compute_mobility(net, mobility, nodeforces, nodetags=[])

   Wrapper to compute nodal velocities.

   :param net: ExaDisNet object
   :param mobility: Mobility object
   :param nodeforces: List of node forces
   :param nodetags: List of NodeTag objects (optional)

.. py:function:: pyexadis.compute_node_mobility(net, i, mobility, fi)

   Wrapper to compute the mobility of a single node.

   :param net: ExaDisNet object
   :param i: Node index (int)
   :param mobility: Mobility object
   :param fi: Force on node

---

Integrator
----------

Integrator Binder
^^^^^^^^^^^^^^^^^

.. py:class:: pyexadis.Integrator

   Integrator wrapper.

   .. py:method:: integrate(system)
      :noindex:

      Integrate the system.

---

Integrator Factories
^^^^^^^^^^^^^^^^^^^^

.. py:function:: pyexadis.make_integrator_trapezoid(params, intparams, force, mobility)

   Instantiate a trapezoid integrator.

.. py:function:: pyexadis.make_integrator_trapezoid_multi(params, intparams, force, mobility)

   Instantiate a multi-step trapezoid integrator.

.. py:function:: pyexadis.make_integrator_rkf(params, intparams, force, mobility)

   Instantiate a RKF integrator.

.. py:function:: pyexadis.make_integrator_rkf_multi(params, intparams, force, mobility)

   Instantiate a multi-step RKF integrator.

.. py:function:: pyexadis.make_integrator_subcycling(params, intparams, force, mobility)

   Instantiate a subcycling integrator.

---

Integration Functions
^^^^^^^^^^^^^^^^^^^^^

.. py:function:: pyexadis.integrate(net, integrator, nodevels, applied_stress, nodetags=[])

   Wrapper to perform a time-integration step.

.. py:function:: pyexadis.integrate_euler(net, params, dt, nodevels, nodetags=[])

   Time-integrate positions using the Euler integrator.

---

Collision
---------

.. py:class:: pyexadis.Collision

   Collision handler.

   .. py:method:: handle(system)
      :noindex:

      Handle collision of the system.

---

.. py:function:: pyexadis.make_collision(collision_mode, params)

   Instantiate a collision class.

.. py:function:: pyexadis.handle_collision(net, collision, xold=[], dt=0.0)

   Wrapper to handle collisions.

---

Topology
--------

.. py:class:: pyexadis.Topology

   Topology handler.

   .. py:method:: handle(system)
      :noindex:

      Handle topology of the system.

---

.. py:function:: pyexadis.make_topology(topology_mode, params, topolparams, force, mobility)

   Instantiate a topology class.

.. py:function:: pyexadis.handle_topology(net, topology, dt)

   Wrapper to handle topological operations.

---

Remesh
------

.. py:class:: pyexadis.Remesh

   Remesh handler.

   .. py:method:: remesh(system)
      :noindex:

      Remesh the system.

---

.. py:function:: pyexadis.make_remesh(remesh_rule, params, remeshparams)

   Instantiate a remesh class.

.. py:function:: pyexadis.remesh(net, remesh)

   Wrapper to remesh the network.

---

Cross-slip
----------

.. py:class:: pyexadis.CrossSlip

   Cross-slip handler.

   .. py:method:: handle(system)
      :noindex:

      Handle cross-slip operations of the system.

---

.. py:function:: pyexadis.make_cross_slip(cross_slip_mode, params, force)

   Instantiate a cross-slip class.

.. py:function:: pyexadis.handle_cross_slip(net, cross_slip)

   Wrapper to handle cross-slip operations.

---


Simulation Driver
-----------------

Drivers
^^^^^^^

.. py:class:: pyexadis.ExaDiSApp

   Base application class for pyexadis simulation driver.

---

.. py:class:: pyexadis.Driver

   Main simulation driver, manages simulation state, modules, and execution.

   .. py:method:: __init__()
      :noindex:

      Default constructor.

   .. py:method:: __init__(system)
      :noindex:

      Construct driver with a system.

      :param system: SystemBind object

   .. py:attribute:: outputdir
      :type: str
      :noindex:

      Output directory path for the simulation.

   .. py:method:: update_state(state)
      :noindex:

      Update the state dictionary with simulation state.

   .. py:method:: read_restart(state, restart_file)
      :noindex:

      Read restart file.

   .. py:method:: set_system(system)
      :noindex:

      Set system for the simulation.

      :param system: SystemBind object

   .. py:method:: set_modules(force, mobility, integrator, collision, topology, remesh, cross_slip=CrossSlipBind())
      :noindex:

      Set modules for the simulation.

      :param force: ForceBind object
      :param mobility: MobilityBind object
      :param integrator: IntegratorBind object
      :param collision: CollisionBind object
      :param topology: TopologyBind object
      :param remesh: RemeshBind object
      :param cross_slip: CrossSlipBind object (optional)

   .. py:method:: set_simulation(restart="")
      :noindex:

      Set things up before running the simulation.

      :param restart: Restart file path (str, optional)

   .. py:method:: initialize(ctrl, check_modules=True)
      :noindex:

      Initialize simulation.

      :param ctrl: Driver.Control object
      :param check_modules: Check modules before initialization (bool, default True)

   .. py:method:: step(ctrl)
      :noindex:

      Execute a simulation step.

   .. py:method:: run(ctrl)
      :noindex:

      Run the simulation.

   .. py:method:: oprec_replay(ctrl, oprec_files)
      :noindex:

      Replay the simulation from OpRec files.

---

Driver Control
^^^^^^^^^^^^^^

.. py:class:: pyexadis.Driver.Control

   Control object for simulation parameters.

   .. py:attribute:: nsteps
      :type: int or Stepper
      :noindex:

      Number of steps or stepper object.

   .. py:attribute:: loading
      :type: Loadings
      :noindex:

      Loading type.

   .. py:attribute:: erate
      :type: float
      :noindex:

      Loading rate.

   .. py:attribute:: edir
      :type: list or array
      :noindex:

      Loading direction.

   .. py:attribute:: appstress
      :type: list or array
      :noindex:

      Applied stress.

   .. py:attribute:: rotation
      :type: bool
      :noindex:

      Enable crystal rotation.

   .. py:attribute:: printfreq
      :type: int
      :noindex:

      Print frequency.

   .. py:attribute:: propfreq
      :type: int
      :noindex:

      Properties output frequency.

   .. py:attribute:: outfreq
      :type: int
      :noindex:

      Configuration and restart output frequency.

   .. py:attribute:: outfreqdt
      :type: float
      :noindex:

      Configuration and restart output time frequency.

   .. py:attribute:: oprecwritefreq
      :type: int
      :noindex:

      OpRec write frequency.

   .. py:attribute:: oprecfilefreq
      :type: int
      :noindex:

      OpRec new file frequency.

   .. py:attribute:: oprecposfreq
      :type: int
      :noindex:

      OpRec nodal positions save frequency.

   .. py:method:: set_props(prop_fields)
      :noindex:

      Set property fields for the output.

---

Driver Loadings Enum
^^^^^^^^^^^^^^^^^^^^

.. py:class:: pyexadis.Driver.Loadings

   Enum for loading types.

   .. py:attribute:: STRESS_CONTROL
      :noindex:

      Stress control loading.

   .. py:attribute:: STRAIN_RATE_CONTROL
      :noindex:

      Strain rate control loading.

---

Driver Stepper
^^^^^^^^^^^^^^

.. py:class:: pyexadis.Driver.Stepper

   Stepper object for simulation steps.

   .. py:method:: __init__(type, stopval)
      :noindex:

      Construct stepper with stopping type and value.

      :param type: Stopping type
      :param stopval: Stopping value

   .. py:method:: NUM_STEPS(nsteps)
      :noindex:

      Construct stepper that iterates to a number of steps.

   .. py:method:: MAX_STEPS(maxsteps)
      :noindex:

      Construct stepper that iterates to a maximum number of steps.

   .. py:method:: MAX_STRAIN(strain)
      :noindex:

      Construct stepper that iterates to a maximum strain.

   .. py:method:: MAX_TIME(time)
      :noindex:

      Construct stepper that iterates to a maximum simulation time.

   .. py:method:: MAX_WALLTIME(time)
      :noindex:

      Construct stepper that iterates to a maximum wall clock time.

   .. py:method:: iterate(driver)
      :noindex:

      Iterate a simulation step.

---
