
ExaDiS Documentation
====================

ExaDiS_ (Exascale Dislocation Simulator) is a set of software modules written to enable numerical simulations of large groups of moving and interacting dislocations, line defects in crystals responsible for crystal plasticity. By tracking time evolution of sufficiently large dislocation ensembles, ExaDiS predicts plasticity response and plastic strength of crystalline materials.

ExaDiS is built around a portable library of core functions for Discrete Dislocation Dynamics (DDD) method specifically written to perform efficient computations on new HPC architectures (e.g. GPUs). Simulations can be driven through the C++ or python interfaces, and interfaced with the OpenDiS_ framework.

.. note:: Although ExaDiS is a fully functional code, it is currently under active development and is subject to frequent updates and bug fixes. There is no guarantee of stability and one should expect occasional breaking changes to the code.

Contributing
------------
As part of the OpenDiS project, ExaDiS is welcoming and encouraging new contributions from users. The easiest way to make contributions to the code is via pull/merge requests, or by reporting bugs or opening discussions in the Issue section of the `ExaDiS Github repo <https://github.com/LLNL/exadis>`_.


License
-------
ExaDiS is released under the BSD-3 license. See LICENSE_ for details.

- Code: LLNL-CODE-862972
- Manual: LLNL-SM-868669


.. _ExaDiS: https://github.com/LLNL/exadis
.. _OpenDiS: https://github.com/OpenDiS/OpenDiS
.. _LICENSE: https://github.com/LLNL/exadis/blob/main/LICENSE


Table of contents
-----------------

.. toctree::
   :maxdepth: 3

   user_guide/index
   developer_guide/index
