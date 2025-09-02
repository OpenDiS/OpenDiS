Python binding: ``pyexadis``
===========================

This section documents the raw ExaDiS classes and functions available through the python interface, ``pyexadis``.
These modules are direct bindings to the backend C++ modules implemented in ExaDiS.
For documentation of the ExaDiS python modules available in OpenDiS format, see :doc:`Python modules: pyexadis_base <python_modules>`.
For documentation about the backend C++ modules please see the :doc:`Developer guide <../../developer_guide/index>` section of the documentation.

.. toctree::
   :maxdepth: 1

Classes
-------

.. py:class:: pyexadis.Params

    ExaDiS simulation parameters.

    .. py:method:: init(crystal="", burgmag, mu, nu, a, maxseg, minseg, rann=-1.0, rtol=-1.0, maxdt=1e-7, nextdt=1e-12, split3node=1) :noindex:

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

    .. py:method:: set_crystal(crystal)

      Set the crystal type.

    .. py:attribute:: crystal :type: str

      Crystal parameters.

    .. py:attribute:: burgmag :type: float

      Burgers vector magnitude (scaling length).

    .. py:attribute:: mu :type: float

      Shear modulus.

    .. py:attribute:: nu :type: float

      Poisson's ratio.

    .. py:attribute:: a :type: float

      Dislocation core radius.

    .. py:attribute:: maxseg :type: float

      Maximum line discretization length.

    .. py:attribute:: minseg :type: float

      Minimum line discretization length.

    .. py:attribute:: rann :type: float

      Annihilation distance.

    .. py:attribute:: rtol :type: float

      Error tolerance.

    .. py:attribute:: maxdt :type: float

      Maximum timestep size.

    .. py:attribute:: nextdt :type: float

      Starting timestep size.

    .. py:attribute:: split3node :type: int

      Enable splitting of 3-nodes.


.. py:class:: pyexadis.CrystalParams

    Parameters for crystal orientation and glide planes.

    .. py:attribute:: R :type: array

      Crystal orientation matrix.

    .. py:attribute:: use_glide_planes :type: bool

      Use and maintain dislocation glide planes.

    .. py:attribute:: enforce_glide_planes :type: bool

      Enforce glide planes option.

    .. py:attribute:: num_bcc_plane_families :type: int

      Number of BCC plane families (1, 2, or 3).


.. py:class:: pyexadis.Crystal

    Crystal type and orientation.

    .. py:method:: init(type) :noindex:

      Initialize with crystal type.

    .. py:method:: init(type, R) :noindex:

      Initialize with crystal type and orientation matrix.

    .. py:attribute:: type :type: int

      Index of the crystal type.

    .. py:attribute:: R :type: array

      Crystal orientation matrix.

    .. py:method:: set_orientation(R)

      Set crystal orientation matrix.

    .. py:method:: set_orientation(euler_angles)

      Set crystal orientation via Euler angles.


.. py:class:: pyexadis.Cell

    Simulation cell and periodic boundary conditions.

    .. py:method:: init(Lbox, centered=False) :noindex:

      Initialize cubic cell.

    .. py:method:: init(Lvecbox, centered=False) :noindex:

      Initialize cell with vector box.

    .. py:method:: init(bmin, bmax) :noindex:

      Initialize cell with bounds.

    .. py:method:: init(h, origin=Vec3(0.0), is_periodic=[PBC_BOUND, PBC_BOUND, PBC_BOUND]) :noindex:

      Initialize cell with matrix, origin, and periodicity.

    .. py:method:: init(cell) :noindex:

      Copy constructor.

    .. py:attribute:: h :type: array

      Cell matrix.

    .. py:attribute:: origin :type: array

      Origin of the cell.

    .. py:method:: center()

      Returns the center of the cell.

    .. py:method:: closest_image(Rref, R)

      Returns the closest image of an array of positions from another reference.

    .. py:method:: pbc_fold(R)

      Fold an array of positions to the primary cell.

    .. py:method:: is_inside(R)

      Checks if a position is inside the primary cell.

    .. py:method:: are_inside(R)

      Checks if an array of positions are inside the primary cell.

    .. py:method:: is_triclinic()

      Returns if the box is triclinic.

    .. py:method:: is_periodic()

      Get the cell periodic boundary condition flags along the 3 dimensions.

    .. py:method:: get_bounds()

      Get the (orthorhombic) bounds of the cell.

    .. py:method:: volume()

      Returns the volume of the cell.


Functions
---------

.. py:function:: pyexadis.initialize(num_threads=-1, device_id=0)

    Initialize the python binding module.

    :param num_threads: Number of OMP threads (default: -1)
    :param device_id: Device (GPU) ID (default: 0)

.. py:function:: pyexadis.finalize()

    Finalize the python binding module.

.. py:function:: pyexadis.read_paradis(filename)

    Read ParaDiS data file.

    :param filename: Path to ParaDiS data file

.. py:function:: pyexadis.generate_prismatic_config(crystal, Lbox, numsources, radius, maxseg=-1, seed=1234, uniform=False)

    Generate a configuration made of prismatic loops in a cubic cell.

    :param crystal: Crystal object
    :param Lbox: Box size
    :param numsources: Number of sources
    :param radius: Loop radius
    :param maxseg: Maximum segment length (default: -1)
    :param seed: Random seed (default: 1234)
    :param uniform: Uniform distribution (default: False)

.. py:function:: pyexadis.generate_prismatic_config(crystal, cell, numsources, radius, maxseg=-1, seed=1234, uniform=False)

    Generate a configuration made of prismatic loops in a custom cell.

    :param crystal: Crystal object
    :param cell: Cell object
    :param numsources: Number of sources
    :param radius: Loop radius
    :param maxseg: Maximum segment length (default: -1)
    :param seed: Random seed (default: 1234)
    :param uniform: Uniform distribution (default: False)
