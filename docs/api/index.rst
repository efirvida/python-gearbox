.. _api:

API reference
=============

Every public class and function of the package, generated from the docstrings. Most programs need
only the models of :mod:`gearbox` (importable from the package itself, ``from gearbox import Gear``)
and the ``calculate_*`` functions of the two standards. The modules whose name starts with an
underscore hold the factors of each standard; they are not part of the stable interface, but they
are listed here because the :doc:`theory pages <../theory/index>` refer to them.

The packages re-export their main names, so the paths of the examples are shorter than the module
paths below:

.. code-block:: python

   from gearbox import Gear, Material, Tool, Transmission         # gearbox.models, gearbox.transmission
   from gearbox.standards import iso6336, agma2101                # iso6336.calculate_pitting, ...
   from gearbox.fem import mesh_gear, solve_bending, solve_contact  # gearbox.fem.mesh, gearbox.fem.calculix

Gear and pair models
--------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   gearbox.models
   gearbox.transmission
   gearbox.enums

Geometry
--------

.. autosummary::
   :toctree: generated
   :nosignatures:

   gearbox.profile
   gearbox.mathutils
   gearbox.optimization

ISO 6336
--------

.. autosummary::
   :toctree: generated
   :nosignatures:

   gearbox.standards.iso6336.pitting
   gearbox.standards.iso6336.bending
   gearbox.standards.iso6336._load_factors
   gearbox.standards.iso6336._materials

AGMA 2001-D04 and AGMA 908-B89
------------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   gearbox.standards.agma2101.pitting
   gearbox.standards.agma2101.bending
   gearbox.standards.agma2101._geometry
   gearbox.standards.agma2101._factors
   gearbox.standards.agma2101._strength
   gearbox.standards._hertz

Finite elements
---------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   gearbox.fem.mesh
   gearbox.fem.calculix
   gearbox.fem.post
