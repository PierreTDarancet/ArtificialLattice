.. TopoQuest documentation master file, created by
   sphinx-quickstart on Tue Jun 25 11:28:06 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TopoQuest's documentation!
=====================================

``TopoQuest`` is a high throughput tight binding based python package the ``quest of topological materials`` with bound states. ``TopoQuest`` is built on ``kwant`` [REF] tight-binding package and ``Z2pack`` package from David Vanderbilt's group. ``TopoQuest`` samples the geometric space of structures satisfying the specified symmetry and geometric constrains as a graph. The graph is passed on to ``kwant`` which constructs the tight binding hamiltonian. The functions in ``Z2pack`` are used to compute the intercell Zak Phase (the origin independant but cell dependant part of Zak Phase) which characterises the topological invariant.    

.. toctree::
   :maxdepth: 2

Code documentation
==================

.. autoclass:: TopoQuest.StructGen.StructGen
    :members:


Redundancy Checker
==================

.. autoclass:: TopoQuest.StructGen.Check_redundant
    :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
