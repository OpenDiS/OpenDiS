PyDiS Documentation
====================

PyDiS (Python Dislocation Simulator) is discrete dislocation dynamics (DDD) simulation module written mostly in Python.  It is developed with two goals in mind.

First, PyDiS provides a ''living'' document about how dislocation dynamics simulation algorithms works.  Both PyDiS and ExaDiS derive from ParaDiS.  There are great similarities in the (high-level) algorithms among these codes, although their implementation details may vary.  So the user can gain a greater understanding of the algorithms by peeking into the implementations of PyDiS, which tend to be pretty short given the expressiveness of Python.

Second, PyDiS offers users an example of rapidly implementing and trying out new algorithms in Python.  Given the modular nature of OpenDiS, the user can run a DDD simulation with some components implemented in Python (e.g. PyDiS) and other components in ExaDiS.  For example, the user may try a new mobility law written in Python while keeping the computationally expensive force calculations in ExaDiS.  
This would greatly shorten the development time while not sacrificing the computational efficiency unnecessarily.  For example, the user may rapidly prototype several mobility laws in Python and select the most desirable one to further implement in ExaDiS style.  By lowering the barrier for prototyping, we expect this set up to greatly shorten the development cycle.

.. note:: PyDiS is currently under active development and is subject to frequent updates and bug fixes. There is no guarantee of stability and one should expect occasional breaking changes to the code.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   simulate_network
   calforce
   mobility_law
   time_integration
   topology
   collision
   remesh




   
