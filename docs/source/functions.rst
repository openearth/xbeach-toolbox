Functions
=========

The XBeach toolbox contains various functions. The functions are defined into the following topics:

General
    General functions, see :ref:`sec_general`

Grid
    Functions related to the grid, see :ref:`sec_grid`

xbeachtools
    Tools to setup a model, see :ref:`sec_xbtools`

xb_analyse
    Tools to analyse a model, see :ref:`sec_XBeachModelAnalysis`




.. _sec_general:
General
--------

Wave Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.dispersion()

.. autofunction:: xbTools.wavecelerity()

.. autofunction:: xbTools.directional_spread_kuik()

.. autofunction:: xbTools.celerity_ratio_equals_09()

.. autofunction:: xbTools.offshore_depth()

Visualize mesh
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.plot_mesh()

.. autofunction:: xbTools.write_mesh_to_shp()

.. autofunction:: xbTools.grid_parameter_to_shape()

.. autofunction:: xbTools.grid_parameter_to_shape()

Executing runs
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.xb_run_script_win()

Geometry
~~~~~~~~~~~~~~~~~~~~


.. autofunction:: xbTools.rotate()

.. autofunction:: xbTools.rotate_grid()

.. autofunction:: xbTools.grid_world2local()

.. autofunction:: xbTools.samples_to_local_grid()

.. autofunction:: xbTools.in_polygon()

.. autofunction:: xbTools.get_points_in_polygon()

.. autofunction:: xbTools.path_distance()



.. _sec_grid:
Grid
---------------------------------------


Refinement
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.grid_refine_grid()

Extionsion
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: xbTools.lateral_extend()

.. autofunction:: xbTools.seaward_extend()


Creation
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.xgrid()

.. autofunction:: xbTools.ygrid()

.. autofunction:: xbTools.grid_transition()



.. _sec_xbtools:
XBeachtools
---------------------------------------

.. autofunction:: xbTools.XBeachModelSetup()



.. _sec_XBeachModelAnalysis:
XBeachModelAnalysis
---------------------------------------

.. autofunction:: xbTools.XBeachModelAnalysis()