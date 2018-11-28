Wraparound: finalize Kwant systems with multiple translational symmetries
=========================================================================

In `Kwant <http://kwant-project.org/>`_ <= 1.2, tight-binding systems with
more than one translational symmetry can be built, but not finalized, due to
limitations of the low-level format.  This module allows to work around this
limitation by replacing all (or all but one) translational symmetries by
parameters (momenta) to the system.  The ensuing ``Builder`` has only 0 or 1 explicit symmetries left and can be thus finalized.  This code makes it possible to

* model transport in systems with periodic boundary conditions of arbitrary
  dimensionality,

* model transport across infinite planes,

* calculate the n-dimensional band-structures.

To test and demo the module, simply execute it: ``python wraparound.py``.  It
works with any Kwant 1.x and both Python 2 and 3.  For instructions, see the
module documentation and the body of the function ``demo()``.

This code was written by Christoph Groth.  The license is the same as that of Kwant (2-clause BSD): http://kwant-project.org/license
