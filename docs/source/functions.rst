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

.. autofunction:: xbTools.general.wave_functions.dispersion()

.. autofunction:: xbTools.general.wave_functions.wavecelerity()

.. autofunction:: xbTools.general.wave_functions.directional_spread_kuik()

.. autofunction:: xbTools.general.wave_functions.celerity_ratio_equals_09()

.. autofunction:: xbTools.general.wave_functions.offshore_depth()

Visualize mesh
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.general.visualize_mesh.plot_mesh()

.. autofunction:: xbTools.general.visualize_mesh.write_mesh_to_shp()

.. autofunction:: xbTools.general.visualize_mesh.grid_parameter_to_shape()

.. autofunction:: xbTools.general.visualize_mesh.grid_parameter_to_shape()

Executing runs
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.general.executing_runs.xb_run_script_win()

Geometry
~~~~~~~~~~~~~~~~~~~~


.. autofunction:: xbTools.general.geometry.rotate()

.. autofunction:: xbTools.general.geometry.rotate_grid()

.. autofunction:: xbTools.general.geometry.grid_world2local()

.. autofunction:: xbTools.general.geometry.samples_to_local_grid()

.. autofunction:: xbTools.general.geometry.in_polygon()

.. autofunction:: xbTools.general.geometry.get_points_in_polygon()

.. autofunction:: xbTools.general.geometry.path_distance()



.. _sec_grid:
Grid
---------------------------------------


Refinement
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.grid.refinement.grid_refine_grid()

Extionsion
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: xbTools.grid.extension.lateral_extend()

.. autofunction:: xbTools.grid.extension.seaward_extend()


Creation
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xbTools.grid.creation.xgrid()

.. autofunction:: xbTools.grid.creation.ygrid()

.. autofunction:: xbTools.grid.creation.grid_transition()



.. _sec_xbtools:
XBeachtools
---------------------------------------

.. autofunction:: xbTools.xbeachtools.XBeachModelSetup()



.. _sec_XBeachModelAnalysis:
XBeachModelAnalysis
---------------------------------------

.. autofunction:: xbTools.xbeachpost.XBeachModelAnalysis()